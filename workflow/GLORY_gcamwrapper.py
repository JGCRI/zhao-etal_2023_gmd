#!/usr/bin/env python
# coding: utf-8

"""
Calculate surface water supply curve for each GCAM basin.

Updated in 2023

Authors:
    Mengqi Zhao (mengqi.zhao@pnnl.gov)
    Thomas Wild (thomas.wild@pnnl.gov)

Project: Global Change Intersectoral Modeling System (GCIMS)

License: BSD 2-Clause
"""

import gcamwrapper
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import copy
from datetime import datetime
from datetime import datetime
from scipy.interpolate import interp1d
from scipy.interpolate import PchipInterpolator
from pyomo.environ import *
from mpi4py import MPI

# =========================================
# |               LP Model                |
# =========================================

def LP_model(K, Smin, Ig, Eg, f, p, z, m):
    """
    Construct Capacity-Yield Curve with Linear Programming Model.

    :param K:       float for reservoir storage capacity
    :param Smin:    int for minimum required active storage for the reservoir
    :param Ig:      float for annual average inflow (over GCAM 5 year time period) to reference reservoirs that has q total storage capacity of Kgref
    :param Eg:      float for annual average evaporation (over GCAM 5 year time period) from reference reservoirs that has q total storage capacity of Kgref
    :param f:       dictionary for demand fraction profile
    :param p:       dictionary for inflow fraction profile
    :param z:       dictionary for evaporation fraction profile
    :param m:       float for fraction of flow released from distributed reservoirs

    :return:        array for capacity - yield curve
    """
    # Linear Programming Model for Water Storage
    model = ConcreteModel()

    # Set model time periods
    TimePeriods = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    LastTimePeriod = 12

    # Define variables
    model.Y = Var(within=NonNegativeReals, initialize=0)
    model.S = Var(TimePeriods, within=NonNegativeReals)
    model.R = Var(TimePeriods, within=NonNegativeReals)
    model.X = Var(TimePeriods, within=NonNegativeReals)
    model.RF = Var(TimePeriods, within=NonNegativeReals)

    # Define parameters
    model.K = Param(initialize=K, within=NonNegativeReals)

    model.f = Param(TimePeriods, initialize=f, within=PercentFraction)
    model.p = Param(TimePeriods, initialize=p, within=PercentFraction)
    model.z = Param(TimePeriods, initialize=z, within=PercentFraction)

    # Define inflow
    def Inflow_init(model, t):
        return max(model.p[t] * (Ig), 0)

    model.I = Param(TimePeriods, initialize=Inflow_init)

    # Define evaporation
    def Evap_init(model, t):
        return max(model.z[t] * (Eg), 0)

    model.E = Param(TimePeriods, rule=Evap_init)

    # Define environmental flow
    def EnvFlow_init(model, t):
        return 0.1 * model.I[t]

    model.EF = Param(TimePeriods, rule=EnvFlow_init)

    # Constraints ----------------------------------------------
    # Define mass balance constraint
    def Mass_rule(model, t):
        if t == LastTimePeriod:
            return model.S[1] == model.S[t] + model.I[t] - model.E[t] - model.EF[t] - model.R[t] + model.RF[t] - \
                   model.X[t]
        else:
            return model.S[t + 1] == model.S[t] + model.I[t] - model.E[t] - model.EF[t] - model.R[t] + model.RF[t] - \
                   model.X[t]

    model.Mass = Constraint(TimePeriods, rule=Mass_rule)

    # Define storage constraint
    def Storage_rule(model, t):
        return (Smin, model.S[t], model.K)

    model.Storage = Constraint(TimePeriods, rule=Storage_rule)

    # Define release constraint
    def Release_rule(model, t):
        return model.R[t] >= model.f[t] * model.Y

    model.Release = Constraint(TimePeriods, rule=Release_rule)

    # Define return flow constraint
    def ReturnFlow_rule(model, t):
        return model.RF[t] == m * (model.R[t] + model.EF[t])

    model.ReturnFlow = Constraint(TimePeriods, rule=ReturnFlow_rule)

    # Objective ----------------------------------------------
    model.obj = Objective(expr=model.Y, sense=maximize)

    # Solve the problem
    opt = SolverFactory('glpk')
    opt.solve(model, tee=False)

    return model.obj()


class DataLoader:
    """
    Load Data
    """

    gcm = 'MIROC-ESM-CHEM'
    rcp = 'rcp6p0'

    data_path = 'your/path/to/inputs_glory'

    f_climate = os.path.join(data_path, 'LP_climate_' + gcm + '_' + rcp + '.csv')
    f_profile = os.path.join(data_path, 'LP_fraction_profile_' + gcm + '_' + rcp + '.csv')
    f_demand = os.path.join(data_path, 'LP_demand_sector_mean_annual_hist.csv')
    f_reservoir = os.path.join(data_path, 'LP_reservoir.csv')
    f_slope = os.path.join(data_path, 'LP_mean_slope_basin.csv')
    f_basin_mapping = os.path.join(data_path, 'basin_to_country_mapping.csv')
    f_region_mapping = os.path.join(data_path, 'basin_to_region_mapping.csv')

    global capacity_gcam_pre

    def __init__(self, basin_id, period, base_period=2020, demand_gcam=None, capacity_gcam=None):
        """
        Initialization

        :param basin_id:            integer for basin id to select
        :param period:              integer for period (year)
        :param base_period:         integer for base year
        :param demand_gcam:         dataframe for GCAM sovled demand
        :param capacity_gcam:       dataframe for storage capacity based on GCAM solved runoff demand
        """

        self.basin_id = basin_id
        self.period = period
        self.base_period = base_period
        self.demand_gcam = demand_gcam
        self.capacity_gcam = capacity_gcam

        self.climate = self.load_data(DataLoader.f_climate)
        self.profile = self.load_data(DataLoader.f_profile)
        self.demand_hist = self.load_data(DataLoader.f_demand)
        self.reservoir = self.load_data(DataLoader.f_reservoir)
        self.slope = self.load_data(DataLoader.f_slope)['slope'].iloc[0]
        self.basin_name_std = self.load_basin_mapping(DataLoader.f_basin_mapping, header_num=7)

        self.inflow = self.climate.loc[self.climate['period'] == self.period, 'smooth_runoff_km3'].iloc[0]
        self.evap_depth = self.climate.loc[self.climate['period'] == self.period, 'pet_km'].iloc[0]
        self.res_area = self.reservoir['nonhydro_area_km2'].iloc[0]

        self.inflow_profile = dict(zip(self.profile.loc[self.profile['period'] == self.period, 'month'],
                                       self.profile.loc[self.profile['period'] == self.period, 'inflow']))
        self.evap_profile = dict(zip(self.profile.loc[self.profile['period'] == self.period, 'month'],
                                     self.profile.loc[self.profile['period'] == self.period, 'pet']))
        self.demand_profile = self.get_demand_profile()

        # capacity values
        self.no_expansion = False
        self.max_capacity = self.get_max_capacity()
        self.current_capacity = self.get_current_capacity()

        # expected incremental size for reservoir storage capacity expansion
        self.expan_incr = self.reservoir['mean_cap_km3'].iloc[0]

        # define constant
        self.storage_min = 0
        self.m = 0.1

    def load_data(self, fn, header_num=0):
        """
        Load data from a CSV file to pandas dataframe.

        :param fn:              string for name of file to load
        :param header_num:      integer for number of lines in file to skip, if text or csv file

        :return:                pandas dataframe
        """
        if not os.path.isfile(fn):
            raise IOError("Error: File does not exist:", fn)

        # for CSV files
        elif fn.endswith('.csv'):
            df = pd.read_csv(fn, skiprows=header_num)
            df = df.loc[df['basin_id'] == self.basin_id]

        else:
            raise RuntimeError("File {} has unrecognized extension".format(fn))

        return df

    def load_gcam_demand(self):
        """
        Load and Format data extract from GCAM using gcamwrapper. The data should be for 235 basins, 6 demand sectors,
        and a single time period.

        :return:                dataframe for sectoral annual demand for selected basin
        """

        if self.demand_gcam is not None:

            df = self.demand_gcam

            # convert objects columns to string
            df = df[['sector', 'subsector', 'year', 'physical-demand']]. \
                rename(columns={'subsector': 'gcam_basin_name', 'physical-demand': 'value'})
            df[['sector', 'gcam_basin_name']] = df[['sector', 'gcam_basin_name']].astype('string', errors='raise')

            # format water withdrawal extracted from GCAM using gcamwrapper
            df['sector'] = df['sector'].str.split(pat='_', expand=True).astype('string', errors='raise').iloc[:, 2]

            # replace demand sector names
            rep = {'an': 'livestock', 'elec': 'electric', 'ind': 'industry',
                   'irr': 'irrigation', 'muni': 'domestic', 'pri': 'mining'}
            df = df.replace({'sector': rep})

            # add standard basin id and basin name
            df = df.merge(self.basin_name_std, how='left', on='gcam_basin_name')

            # filter to selected basin
            df = df.loc[df['basin_id'] == self.basin_id]

            # aggregate demand for each sector within each basin
            grp = df.groupby(['basin_id', 'gcam_basin_name', 'sector', 'year'], as_index=False).sum()

            # rename column 'value'
            grp = grp.rename(columns={'value': 'demand_ann'})

        else:
            grp = None

        return grp

    @staticmethod
    def load_basin_mapping(fn, header_num=7):
        """
        Mapping different formats of basin names.

        :param fn:              string for full file path
        :param header_num:      integer for numbers of rows to skip until the header

        :return:                dataframe
        """

        # load basin mapping data
        df_basin = pd.read_csv(DataLoader.f_basin_mapping, skiprows=header_num)
        df_region = pd.read_csv(DataLoader.f_region_mapping)

        # select relevant columns and rename
        df_basin = df_basin.loc[:, ['GCAM_basin_ID', 'Basin_long_name', 'GLU_name']]

        df_basin = df_basin.rename(columns={'GCAM_basin_ID': 'basin_id',
                                            'Basin_long_name': 'basin_name',
                                            'GLU_name': 'gcam_basin_name'})

        # convert from object to string or int
        df_basin[['basin_name', 'gcam_basin_name']] = df_basin[['basin_name', 'gcam_basin_name']].astype('string', errors='raise')
        df_basin['basin_id'] = df_basin['basin_id'].astype('int32', errors='raise')

        # join region name
        df = pd.merge(df_basin, df_region, how='left', on=['gcam_basin_name'])

        return df

    def get_demand_profile(self):
        """
        Calculate total demand profile with historical sectoral profile and sectoral demand.

        :param data_type:           string for input demand data. 'hist' or 'gcam'

        :return:                    dictionary
        """

        if self.period <= self.base_period:
            # get historical demand
            demand_ann = self.demand_hist[['sector', 'demand_ann']]

        elif self.period > self.base_period:
            if self.load_gcam_demand() is None:
                demand_ann = self.demand_hist[['sector', 'demand_ann']]
            else:
                if self.load_gcam_demand()['demand_ann'].sum() == 0:
                    print('Basin: ', self.basin_id, ' has a sum of 0 demand from all sectors. Replace demand profile with historical profile.')
                    demand_ann = self.demand_hist[['sector', 'demand_ann']]
                else:
                    # reformat gcam withdrawal
                    demand_ann = self.load_gcam_demand()[['sector', 'demand_ann']]

        # only keep sectoral demand profiles and melt
        df = self.profile.loc[self.profile['period'] == self.period].copy()
        df = df.drop(['basin_id', 'basin_name', 'period', 'pet', 'inflow'], axis=1). \
            melt(id_vars=['month']).rename(columns={'variable': 'sector'})

        # merge annual demand and sectoral demand profiles
        df = pd.merge(df, demand_ann, how='left', on=['sector'])

        # calculate demand amount for each demand sector
        df['demand'] = df['value'] * df['demand_ann']

        # calculate total monthly demand by aggregating all demand sectors
        df = df.groupby('month', as_index=False).sum()

        # calculate profile
        df['profile'] = df['demand'] / df['demand_ann']

        # construct dictionary for month and demand profile
        dict_out = dict(zip(df['month'], df['profile']))

        return dict_out

    def get_current_capacity(self):
        """
        calculate if there will be expansion on storage capacity.

        :return:                float value for storage capacity
        """

        if self.period <= self.base_period:
            # existing storage capacity
            val = self.reservoir['nonhydro_cap_km3'].iloc[0]

        if self.capacity_gcam is not None:
            if self.period == self.base_period + 5:

                val = max(self.reservoir['nonhydro_cap_km3'].iloc[0],
                          self.capacity_gcam.loc[self.capacity_gcam['basin_id'] == self.basin_id, 'capacity_gcam'].values[0])
            else:
                # if gcam solved capacity is larger than previous capacity, then expand
                capacity_gcam_pre_basin = capacity_gcam_pre.loc[capacity_gcam_pre['basin_id'] == self.basin_id, 'capacity_gcam'].values[0]
                capacity_gcam_basin = self.capacity_gcam.loc[self.capacity_gcam['basin_id'] == self.basin_id, 'capacity_gcam'].values[0]
                if capacity_gcam_basin > capacity_gcam_pre_basin:
                    val = capacity_gcam_basin
                else:
                    val = capacity_gcam_pre_basin
        else:
            val = self.reservoir['nonhydro_cap_km3'].iloc[0]

        return val

    def get_max_capacity(self):
        """
        Adjust maximum storage capacity value only if the basin have no expandable capacity.
        Max storage capacity input data for some basins is already adjusted based on the historical storage capacity.

        :return:            float64
        """
        val = self.reservoir['expan_cap_km3'].iloc[0]

        # adjust max storage cap to mean annual runoff both max and current capacities are 0
        if val == 0:
            val = 0.01 * self.climate.loc[self.climate['period'] == self.period, 'smooth_runoff_km3'].iloc[0]

        return val


class SupplyCurve:
    """
    Calculate supply curves based on Capacity-Yield curves from linear programming model

    USAGE:      SupplyCurve(basin_id=id,
                            period=period,
                            demand_gcam=demand_gcam,
                            capacity_gcam=capacity_gcam)

    """

    global fig_path

    def __init__(self, basin_id, period, demand_gcam, capacity_gcam):
        self.basin_id = basin_id
        self.period = period
        self.demand_gcam = demand_gcam
        self.capacity_gcam = capacity_gcam

        self.d = DataLoader(basin_id=self.basin_id,
                            period=self.period,
                            demand_gcam=self.demand_gcam,
                            capacity_gcam=self.capacity_gcam)

        # initialized base yield when capacity is 0
        self.yield_base = self.run_lp_model(capacity=0, capacity_track=0)

        # initialize output list
        self.output = []
        self.output.append([0, self.yield_base])

        # initialize
        self.init_intv = 100  # initial number of intervals for reservoir storage capacity expansion
        self.discount_rate = 0.05  # discount rate for reservoir capital cost
        self.life_time = 60  # lifetime in years
        self.base_price = 1e-04  # starting price for water, $/m3
        self.mantainance = 0.1  # maintanance fee as fraction of reservoir cost

        # parameter
        self.EXCEED_ADJ = 0.8  # multiplier that adjusts variables that goes over maximum allowable values
        self.EXCEED_ADJ_ACTION = False  # binary that shows whether there is an adjustment for separation point for expansion sequence

        self.dx = None  # increment for reservoir storage capacity expansion
        self.inflection_point = None  # inflection point on capacity-yield curve
        self.max_expansion_point = None  # the capacity and yield values for the max expansion point
        self.expansion_unit = None  # storage capacity increments
        self.max_supply = None  # maximum supply
        self.cost_per_m3 = None # reservoir capital cost in 1975USD/m3

        # discrete capacity-yield curve unconstrained by maximum capacity
        self.capacity_yield_unconstrained = self.get_capacity_yield_unconstrained()

        # interpolatable capacity-yield curve
        self.capacity_yield_interpld = interp1d(self.capacity_yield_unconstrained['capacity'],
                                                self.capacity_yield_unconstrained['yield'],
                                                bounds_error=False,
                                                fill_value=max(self.capacity_yield_unconstrained['yield']))

        # get capacity yield curve (constrained by maximum capacity)
        self.capacity_yield = self.get_capacity_yield_constrained()

        # storage capacity expansion sequence
        self.expansion_seq = self.get_expansion_sequence()

        # yield increment sequence based on storage capacity expansion path
        self.yield_increment_seq = self.get_yield_gain_sequence()

        # cost per unit expansion
        self.cost = self.cost_per_expansion()

        # price
        self.price = self.levelized_cost()

        # supply curve
        self.supply_curve = self.construct_supply_curve()

        # maxsubresource
        self.maxsubresource = self.construct_maxsubresource()

        # plot
        self.plot_curve(save_fig=True)

    def get_reservoir_evap(self, capacity_track):
        """
        calculate reservoir surface area based on current storage capacity.
        Calculate reservoir ET in km3/year based on the new storage capacity in each basin.

        :return:
        """

        if isinstance(capacity_track, float) or isinstance(capacity_track, int):
            df_ini = pd.DataFrame({capacity_track}, columns=['capacity'])
        elif isinstance(capacity_track, list):
            df_ini = pd.DataFrame(capacity_track, columns=['capacity'])
        elif capacity_track is None:
            raise TypeError('Capacity Track can not be None.')

        start = df_ini.iloc[0, 0]
        end = df_ini.iloc[-1, 0]

        if end - start > self.d.expan_incr:
            df = pd.DataFrame(np.append(np.arange(start, end, self.d.expan_incr), end), columns=['capacity'])
        else:
            df = df_ini.copy()

        # calculate the changes in capacity
        df['increase'] = df['capacity'].diff()
        df.loc[0, 'increase'] = df.loc[0, 'capacity']

        # get b and c for Volume-Area Correlation V = c * A^b
        b = self.d.reservoir['b'].iloc[0]
        c = self.d.reservoir['c'].iloc[0]

        # calculate corresponding area for each increase in capacity
        df['area'] = (df['increase'] / c) ** (1 / b)

        # calculate et volume
        df['evap'] = df['area'] * self.d.evap_depth

        # calculate total evaporation from total capacity
        evap_vol = sum(df['evap'])

        if evap_vol > self.d.inflow:
            evap_vol = self.d.inflow * self.EXCEED_ADJ

        return evap_vol

    def run_lp_model(self, capacity, capacity_track):
        """
        Run LP model.

        :param capacity:                float for target capacity
        :param capacity_track           dataframe for tracking capacity at each iteration

        :return:                        value for optimized yield
        """

        evap_vol = self.get_reservoir_evap(capacity_track=capacity_track)

        val = LP_model(K=capacity, Smin=self.d.storage_min, Ig=self.d.inflow, Eg=evap_vol,
                       f=self.d.demand_profile, p=self.d.inflow_profile, z=self.d.evap_profile, m=self.d.m)

        return val

    def set_dx(self):
        """
        The algorithm of choosing dx in this section is to prevent the LP model running over 500 iterations,
        or not enough iterations (<50).
        However, another thing to notice is that the inflection point when reaches mean annual inflow.
        This requires we make sure there are enough iteration before the inflection point (20 or more).
        The algorithm is used to detect if the yield value with K = dx is less than 0.05 of Ig or
        Y0+0.05(Ig-Y0) if Y0 != 0, dx is reduce to 1/2 until dx satisfies the condition.

        :return:                none
        """
        
        dx_coarse = self.d.max_capacity / self.init_intv
        iteration = int(math.ceil(self.d.max_capacity / self.d.expan_incr))

        if (iteration > 500) or (iteration < 50):
            self.dx = dx_coarse
        else:
            self.dx = self.d.expan_incr

        # initialize init yield from first incremental capacity expansion
        yield_dx = self.run_lp_model(capacity=self.dx, capacity_track=self.dx)

        # deduce dx so that it is small enough where the increase in yield value is less than 5% of the difference
        # between average annual inflow and base yield
        n = 0
        while yield_dx > (self.yield_base + 0.05 * (self.d.inflow - self.yield_base)):
            self.dx = self.dx / 2
            yield_dx = self.run_lp_model(capacity=self.dx, capacity_track=self.dx)
            n += 1

        # increase dx so that the initial incremental yield from base yield if capacity expand dx is large enough
        # to make half of the intervals before inflection point
        n = 0
        try:
            while (self.d.inflow - self.yield_base) / (yield_dx - self.yield_base) > self.init_intv / 2:
                #print('Calculating Basin - ', self.basin_id, ', current dx is ', self.dx)
                self.dx = self.dx * 2
                yield_dx = self.run_lp_model(capacity=self.dx, capacity_track=self.dx)
                n += 1
        except ZeroDivisionError as err:
            pass

    def get_capacity_yield_unconstrained(self):
        """
        Calculate capacity yield curve.

        :return:                dataframe
        """
        # set dx for the curve before inflection point
        self.set_dx()

        # calculate capacity-yield curve before reaching inflection point
        i = 0
        capacity_cum = 0
        capacity_track = [capacity_cum]
        while self.output[i][1] < (self.yield_base + 0.95 * (self.d.inflow - self.yield_base)):
            capacity_cum += self.dx
            capacity_track.append(capacity_cum)
            yield_cum = self.run_lp_model(capacity=capacity_cum, capacity_track=capacity_track)

            if yield_cum - self.output[i][1] > 0:
                self.output.append([capacity_cum, yield_cum])
                i += 1
            else:
                capacity_track.remove(capacity_track[-1])
                capacity_cum = capacity_track[-1]
                break

        # average slope from the last 5 points from the inflection point
        slope = (np.diff(list(list(zip(*self.output[-5:]))[1])) / np.diff(list(list(zip(*self.output[-5:]))[0]))).mean()

        # determine the dx for the part approaching inflection point to be enough to capture the curve details
        while self.output[i][1] - self.output[i - 1][1] > 0:
            if slope < 1.5:
                capacity_cum += min(self.dx / 2, self.d.expan_incr)
                capacity_track.append(capacity_cum)
            else:
                capacity_cum += min(self.dx / 10, self.d.expan_incr)
                capacity_track.append(capacity_cum)

            yield_cum = self.run_lp_model(capacity=capacity_cum, capacity_track=capacity_track)

            if yield_cum > self.output[i][1]:
                self.output.append([capacity_cum, yield_cum])
                i += 1
            else:
                capacity_track.remove(capacity_track[-1])
                capacity_cum = capacity_track[-1]
                print('Basin', self.basin_id, ':', i, 'iterations')
                break

        # inflection point is the point on the capacity-yield curve where yield will NOT increase with capacity
        self.inflection_point = pd.DataFrame({'capacity': [self.output[-1][0]], 'yield': [self.output[-1][1]]})

        if self.output[-1][0] < self.d.max_capacity:
            # if maximum expandable capacity is higher than the inflection point
            # calculate the yield at max expandable capacity
            # reservoir construction after inflection point is still by capacity_cum size
            capacity_extend = np.append(
                np.arange(self.output[-1][0] + capacity_cum, self.d.max_capacity, capacity_cum),
                self.d.max_capacity)
            capacity_track.extend(capacity_extend.tolist())
            yield_cum = self.run_lp_model(capacity=self.d.max_capacity, capacity_track=capacity_track)

            self.max_expansion_point = pd.DataFrame({'capacity': [self.d.max_capacity],
                                                     'yield': [yield_cum],
                                                     'exceed_inflection': 'yes'})

        return pd.DataFrame(self.output, columns=['capacity', 'yield'])

    def get_capacity_yield_constrained(self):
        """
        Constrain the capacity yield curve by max capacity.

        :return:
        """

        # # if the inflection point goes beyond max capacity, keep the curve to max capacity
        # capacity_yield_constrained = self.capacity_yield_unconstrained.loc[
        #     self.capacity_yield_unconstrained['capacity'] < self.d.max_capacity]
        #
        # # if the inflection point goes beyond max capacity, keep the curve after max capacity
        # capacity_yield_exceed = self.capacity_yield_unconstrained.loc[
        #     self.capacity_yield_unconstrained['capacity'] > self.d.max_capacity]

        # get maximum yield constrained by maximum capacity
        # This value may be adjusted if max cap > inflection point
        max_expansion_yield_adj = self.capacity_yield_interpld(self.d.max_capacity)[()]

        max_cap_yield = pd.DataFrame({'capacity': [self.d.max_capacity], 'yield': [max_expansion_yield_adj]})

        if self.max_expansion_point is None:
            self.max_expansion_point = pd.DataFrame({'capacity': [self.d.max_capacity],
                                                     'yield': [max_expansion_yield_adj],
                                                     'exceed_inflection': 'no'})

        # append constrained max point
        capacity_yield_constrained = pd.concat([self.capacity_yield_unconstrained, max_cap_yield], ignore_index=True).\
            sort_values(by=['capacity', 'yield'], ignore_index=True).drop_duplicates()

        capacity_yield_constrained.reset_index(drop=True, inplace=True)

        return capacity_yield_constrained

    def get_expansion_sequence(self):
        """
        Determine the evenly spaced intervals for incremental storage capacity and get capacity expansion sequence.

        :return:                1D array
        """
        # find the first index when it just reaches the highest yield
        index_inflection = self.capacity_yield[self.capacity_yield['yield'] == max(self.capacity_yield.iloc[:, 1])].index[0]

        # Method 1 ============
        # # determine numbers of intervals
        # num = max(5, int(math.ceil(self.capacity_yield.iloc[index_inflection, 0] / self.d.expan_incr)))
        #
        # # capacity expansion increments
        # self.expansion_unit = self.capacity_yield.iloc[index_inflection, 0] / num
        #
        # # set the sequence of storage capacity expansion
        # expansion = np.linspace(0, self.capacity_yield.iloc[index_inflection, 0], num)

        # Method 2 ============
        # set expansion unit the same as defined expan_incr
        self.expansion_unit = self.d.expan_incr

        # capacity at the inflection point
        cap_inflection = self.capacity_yield.iloc[index_inflection, 0]

        # initial array with interval of expan_incr
        expansion = np.arange(0,
                              cap_inflection,
                              self.expansion_unit)

        # append the last point if the value is not extremely small
        if (cap_inflection - expansion[-1]) / self.expansion_unit > 0.5:
            expansion = np.append(expansion, cap_inflection)

        return expansion

    def get_yield_gain_sequence(self):
        """
        Calculate yield increments sequence based on reservoir storage capacity expansion.
        two segments of yield increment sequences: (1) before max capacity and (2) after max capacity

        :return:                1D array
        """

        # determine the break point to separate expansion sequence
        # this is to avoid the situation when max capacity exceed inflection point, because
        # 1. it does not make sense to keep expanding capacity when no more water can be provided
        # 2. the expansion sequence after separation point means using other facilities rather than reservoirs,
        #    and the costs will be higher to create enough curvature for supply curves
        if self.inflection_point.loc[0, 'capacity'] < self.d.max_capacity:
            separate_point = self.inflection_point.loc[0, 'capacity'] * self.EXCEED_ADJ
            self.EXCEED_ADJ_ACTION = True
        else:
            separate_point = self.d.max_capacity

        expansion_seq_pre = self.expansion_seq[self.expansion_seq <= separate_point]
        expansion_seq_after = self.expansion_seq[self.expansion_seq > separate_point]

        yield_gain_pre = self.capacity_yield_interpld(expansion_seq_pre[1:]) - \
                         self.capacity_yield_interpld(expansion_seq_pre[:-1])

        if len(expansion_seq_after) == 0:
            yield_gain_after = None
        else:
            expansion_seq_after = np.insert(expansion_seq_after, 0, expansion_seq_pre[-1])
            yield_gain_after = self.capacity_yield_interpld(expansion_seq_after[1:]) - \
                               self.capacity_yield_interpld(expansion_seq_after[:-1])

        yield_gain = [yield_gain_pre, yield_gain_after]

        return yield_gain

    def cost_per_expansion(self):
        """
        Calculate the cost in 1975 USD per unit reservoir storage capacity expansion

        :return:                float in million 1975$
        """
        # cost = 15.375 * pow(Kexp, 0.3995) # million AUD$ For watersupply reservoirs, unit GL (Petheram et al., 2019)
        # convert from GL to unit km3 (Petheram et al., 2019)
        # convert from AUD to USD for 2016: 1 USD = 1.345 AUD in 2016
        # (source: https://data.worldbank.org/indicator/PA.NUS.FCRF?end=2020&locations=AU&start=2006)
        # convert 2016 USD to 1975 USD = 101.049/28.511 = 3.544: Deflators 1975 (28.511) to 2016 (101.049)
        # (source: World Bank https://data.worldbank.org/indicator/NY.GDP.DEFL.ZS?locations=US&view=chart)

        # val = 242.8371 * pow(self.expansion_unit, 0.3995) / 1.345 / 3.544

        # convert capacity from km3 to million m3
        cap_mcm = self.expansion_unit * 1e3

        # slope
        x = self.d.slope

        # calculate estimated cost based on storage cose - slope relationships
        # this is normalized cost in cost/m3, the value still need to be converted to 1975 USD
        if 0 <= cap_mcm < 25:
            val = 0.0197 * x ** 2 + 0.0538 * x + 0.5818
        elif 25 <= cap_mcm < 49:
            val = 0.0295 * x ** 2 - 0.0044 * x + 0.4456
        elif 49 <= cap_mcm < 74:
            val = 0.034 * x ** 2 - 0.031 * x + 0.3982
        elif 74 <= cap_mcm < 123:
            val = 0.037 * x ** 2 - 0.0521 * x + 0.3655
        elif 123 <= cap_mcm < 247:
            val = 0.0372 * x ** 2 - 0.0607 * x + 0.3094
        elif 247 <= cap_mcm < 493:
            val = 0.0368 * x ** 2 - 0.0671 * x + 0.2633
        elif 493 <= cap_mcm < 1233:
            val = 0.0372 * x ** 2 - 0.0607 * x + 0.3094
        elif 1233 <= cap_mcm < 2467:
            val = 0.0362 * x ** 2 - 0.0824 * x + 0.1895
        elif 2467 <= cap_mcm < 4934:
            val = 0.0368 * x ** 2 - 0.0671 * x + 0.2633
        elif 4934 <= cap_mcm < 12335:
            val = 0.0334 * x ** 2 - 0.0868 * x + 0.1427
        elif 12335 <= cap_mcm:
            val = 0.0314 * x ** 2 - 0.0896 * x + 0.1111

        # assuming the average unit cost of reservoirs based on
        # Keller, A., R. Sakthivadivel and D. Seckler. 2000. “Water Scarcity and the Role of Storage in Development.”
        # estimate mean of small and medium project (~ 0.13 * 4 = 0.52 1998$/m3 ),
        # and mean of large storage projects (~ 0.11 * 4 = 0.44 1998$/m3) and take mean between two
        # Table 5 is in 1998 US dollars and need to be converted to 1975 USD
        # https://data.worldbank.org/indicator/FP.CPI.TOTL?locations=US
        mean_global_reservoir_cost = (0.52 + 0.44) / 2

        convert_1998_to_1975_usd = 24.7 / 74.8

        # final unit is million 1975$ per expansion unit
        cost_per_expansion = (1 + self.mantainance) * val * self.expansion_unit * 10**9 * mean_global_reservoir_cost * convert_1998_to_1975_usd / 10**6

        # calculate cost per m3
        self.cost_per_m3 = np.array([[self.basin_id, cost_per_expansion * 10e6 / (self.expansion_unit * 10**9)]])

        return cost_per_expansion

    def levelized_cost(self):
        """
        Calculate levelized cost of storage capacity (LCOSC) as 1975 USD per unit water yield.

        :return:                1D numpy array
        """
        if self.EXCEED_ADJ_ACTION:
            cost_arr = np.array([self.cost, self.cost * 10])
        else:
            cost_arr = np.array([self.cost, self.cost * 5])

        with np.errstate(divide='ignore'):
            price_increment_pre = np.where(self.yield_increment_seq[0] == 0,
                                           0,
                                           (cost_arr[0] * 10 ** 6 * self.discount_rate /
                                            (1 - pow(1 + self.discount_rate, -self.life_time))) / (
                                                   self.yield_increment_seq[0] * 10 ** 9))
            if self.yield_increment_seq[1] is not None:
                price_increment_after = np.where(self.yield_increment_seq[1] == 0,
                                                 0,
                                                 (cost_arr[1] * 10 ** 6 * self.discount_rate /
                                                  (1 - pow(1 + self.discount_rate, -self.life_time))) / (
                                                         self.yield_increment_seq[1] * 10 ** 9))
                price_increment = np.concatenate((price_increment_pre, price_increment_after), axis=None)
            else:
                price_increment = price_increment_pre

            # when there is no reservoir, water costs $0
            price_increment = np.insert(price_increment, 0, self.base_price)

            price = np.cumsum(price_increment, dtype=float)

        return price

    def get_20_point_supply_curve(self, supply_curve_raw):
        """
        This is to get 20-point sequence of yields with monotonic cubic spline.

        :param supply_curve_raw:            dataframe for raw supply curve without supply sequence spacing and smoothing

        :return:
        """
        # monotonic cubic spline interpolation
        supply_curve_interpld = PchipInterpolator(supply_curve_raw['supply'], supply_curve_raw['price'])

        seq_sample = np.linspace(0, self.max_supply, num=100, endpoint=True)

        # calculate 2nd order derivative
        dx = seq_sample[1] - seq_sample[0]
        y = supply_curve_interpld(seq_sample)
        dydx = np.gradient(y, dx)

        # mean of the gradient
        grad_mean = np.mean(dydx)

        # inflection index
        index_inflection = np.where(dydx < grad_mean)[0][-1]

        # keep the high gradient part of the sequence
        index_back = np.where(dydx > grad_mean)[0]
        index_back = index_back[index_back > index_inflection]
        if len(index_back) or (index_back[-1] - index_back[0]) > 10:
            space = max(math.ceil(len(index_back) / 10), math.ceil((index_back[-1] - index_back[0]) / 10))
            index_back = np.unique(np.append(np.arange(index_back[0], index_back[-1] + 1, space), index_back[-1]))

        # the first two front indexes need to be at lease 0 and 1e-5, to have enough curvature for GCAM
        price_seq_front = supply_curve_interpld(seq_sample[0:index_inflection])
        price_index = np.where(price_seq_front > self.base_price)[0]

        if len(price_index) == 0:
            price_index = np.where(price_seq_front > self.base_price / 10)[0]

        # calculate space of the flat part of the sequence
        space = max(math.ceil((price_index[-1] - price_index[0])/(17 - len(index_back))), 1)

        # get flat part sequence index
        index_front = np.arange(price_index[0], price_index[-1] + 1, space)
        index_front = np.insert(index_front, 0, 0)

        # check if total length of front and back index is smaller than 18
        i = 1
        length = len(index_front) + len(index_back)
        while len(index_front) + len(index_back) < 18:
            space = math.ceil((price_index[0] / 3) / (18 - length))
            index_front = np.append(index_front, price_index[0] - space * i)
            i += 1

        # full index
        index = np.sort(np.concatenate((index_front, index_back), axis=None))

        # sequence
        seq = seq_sample[index]

        # interpolate price in supply curve
        supply_curve = pd.DataFrame({'supply': seq, 'price': supply_curve_interpld(seq)})

        # smooth the curve with rolling mean, which makes curve 20 point in most cases
        window = 5
        if window > 0:
            supply_curve_forward = supply_curve.rolling(window=window, min_periods=1).mean()
            supply_curve_backward = supply_curve[::-1].rolling(window=window, min_periods=1).mean()[::-1]
            supply_curve = pd.concat([supply_curve_forward, supply_curve_backward]).round(decimals=6)
            supply_curve = supply_curve.drop_duplicates().sort_values(by=['supply', 'price'], ignore_index=True)
            
        supply_curve_interpld = PchipInterpolator(supply_curve['supply'], supply_curve['price'])

        j = 2
        while len(supply_curve) > 20:
            supply_curve.drop(supply_curve.index[-j], axis=0, inplace=True)
            j += 2

        while len(supply_curve) < 20:
            seq_append = seq_sample[np.where(dydx < grad_mean)[0][-i]]
            temp = pd.DataFrame({'supply': [seq_append], 'price': [supply_curve_interpld(seq_append)[()]]})
            supply_curve = pd.concat([supply_curve, temp]).sort_values(by=['supply', 'price'], ignore_index=True)
            i += 2

        supply_curve.reset_index(drop=True, inplace=True)

        return supply_curve

    def get_20_even_point_supply_curve(self, supply_curve_raw):
        """
        This is to get 20 evenly spaced point on supply sequence to get the curve.

        :param supply_curve_raw:
        :return:
        """

        supply_curve = supply_curve_raw

        # get the 20 evenly spaced points for x axis (supply)
        supply_seq = np.linspace(supply_curve.loc[0, 'supply'], self.max_supply, num=20, endpoint=True)

        # smooth the curve with rolling mean, which makes curve 20 point in most cases
        window = 5
        if window > 0:
            supply_curve_forward = supply_curve.rolling(window=window, min_periods=1).mean()
            supply_curve_backward = supply_curve[::-1].rolling(window=window, min_periods=1).mean()[::-1]
            supply_curve = pd.concat([supply_curve_forward, supply_curve_backward]).round(decimals=6)
            supply_curve = supply_curve.drop_duplicates().sort_values(by=['supply', 'price'], ignore_index=True)

        # Interpolate smoothed supply curve
        supply_curve_interpld = PchipInterpolator(supply_curve['supply'], supply_curve['price'])

        # 20 points supply curve
        supply_curve_20pt = pd.DataFrame({'supply': supply_seq, 'price': supply_curve_interpld(supply_seq)})

        return supply_curve_20pt


    def construct_supply_curve(self):
        """
        Construct supply curve.
        
        :return:                2D numpy array
        """
        # supply for all expansion points
        supply = self.capacity_yield_interpld(
            np.unique(np.concatenate((self.expansion_seq[:-1], self.expansion_seq[1:]))))

        # get maximum supply
        self.max_supply = supply[-1]

        # combine supply and price in to 2D array
        supply_curve_lp = np.vstack((supply, self.price)).T

        # starting point for supply curve
        supply_curve_min = np.array([0, 0], ndmin=2)

        # defined the entire curve
        supply_curve = pd.DataFrame(np.concatenate((supply_curve_min, supply_curve_lp), axis=0),
                                    columns=['supply', 'price'])

        # 20 points supply curve
        supply_curve_20pt = self.get_20_point_supply_curve(supply_curve_raw=supply_curve)
        # supply_curve_20pt = self.get_20_even_point_supply_curve(supply_curve_raw=supply_curve)

        # update max supply since the digit is limited to 6
        self.max_supply = supply_curve_20pt.iloc[19, 0]

        # Append 20 points supply curve to data frame
        gcam_basin_name = \
        self.d.basin_name_std.loc[self.d.basin_name_std['basin_id'] == self.basin_id, ['gcam_basin_name']].iloc[0, 0]
        supply_curve_append = pd.DataFrame({'resource': gcam_basin_name + '_water withdrawals',
                                            'grade': ['grade{}'.format(i) for i in range(1, 21)],
                                            'available': supply_curve_20pt['supply'] / self.max_supply,
                                            'extractioncost': supply_curve_20pt['price']})

        # change to numpy
        supply_curve_array = supply_curve_append.to_numpy()

        return supply_curve_array

    def construct_maxsubresource(self):
        """
        Construct new maxSubResource dataframe based on the max supply constrained by the max capacity.

        :return:
        """

        # extract region that has the largest portion of the basin
        region = self.d.basin_name_std.loc[self.d.basin_name_std['basin_id'] == self.basin_id, 'region'].iloc[0]

        # renewresource
        resource = self.d.basin_name_std.loc[
            self.d.basin_name_std['basin_id'] == self.basin_id, 'gcam_basin_name'].iloc[0] + '_water withdrawals'

        maxsubresource = pd.DataFrame({'region': [region],
                                       'resource': [resource],
                                       'subresource': ['runoff'],
                                       'year': [self.period],
                                       'maxSubResource': [self.max_supply]})

        # change to numpy
        maxsubresource_array = maxsubresource.to_numpy()

        return maxsubresource_array

    def plot_curve(self, save_fig=True):
        """
        Plot capacity yield curve.
        """
        basin_name = self.d.basin_name_std.loc[self.d.basin_name_std['basin_id'] == self.basin_id, 'basin_name'].iloc[0]

        fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 8), constrained_layout=False)

        # plot capacity - yield curve
        ax[0].plot(self.capacity_yield['capacity'], self.capacity_yield['yield'], color='black', lw=3)
        max_cap, = ax[0].plot(self.d.max_capacity, self.capacity_yield_interpld(self.d.max_capacity),
                              linestyle='None', marker='o', markersize=14, markeredgecolor='black',
                              markerfacecolor='darkorange')
        current_cap, = ax[0].plot(self.d.current_capacity, self.capacity_yield_interpld(self.d.current_capacity),
                                  linestyle='None', marker='o', markersize=10, markeredgecolor='black',
                                  markerfacecolor='lightseagreen')
        ax[0].set_title(str(self.basin_id) + ' - ' + basin_name, fontsize=16, fontweight='bold')
        ax[0].legend([max_cap, current_cap], ['max expandable capacity', 'current capacity'], numpoints=1)
        ax[0].set_xlabel('Storage Capacity ($km^3$)', fontsize=14)
        ax[0].set_ylabel('Yield ($km^3/year$)', fontsize=14)
        ax[0].grid()

        # plot supply curve
        ax[1].plot(self.supply_curve[:, 2], self.supply_curve[:, 3], color='black', lw=3,
                   marker='o', markersize=8, markeredgecolor='black', markerfacecolor='dodgerblue')
        # ax[1].set_title(str(self.basin_id) + ' - ' + basin_name, fontsize=16, fontweight='bold')
        ax[1].set_xlabel('Availability', fontsize=14)
        ax[1].set_ylabel('Price (1975USD/year)', fontsize=14)
        ax[1].grid()

        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        fig.tight_layout()

        if save_fig:
            fig_path_period = os.path.join(fig_path, str(self.period))
            fig.savefig(os.path.join(fig_path_period, str(self.basin_id) + ' - ' + basin_name + '.png'))

        # plt.show()


class InfoExchange:
    """
    Pass infor between GCAM and LP model.
    """
    SUPPLY_COST_QUERY = "world/region{region@name}/resource[+NamedFilter,StringRegexMatches,_water withdrawals]/subresource[+NamedFilter,StringEquals,runoff]/grade{grade@name}/extractioncost"
    SUPPLY_AVAIL_QUERY = "world/region{region@name}/resource[+NamedFilter,StringRegexMatches,_water withdrawals]/subresource[+NamedFilter,StringEquals,runoff]/grade{grade@name}/available"
    MAXSUBRESOURCE_QUERY = "world/region{region@name}/resource[+NamedFilter,StringRegexMatches,_water withdrawals]/subresource[+NamedFilter,StringEquals,runoff]/maxSubResource{year@year}"
    PRICE_SOLVE_QUERY = "marketplace/market[+NamedFilter,StringRegexMatches,_water withdrawals]/market-period{year@year}/price"
    WITHDRAWAL_SECTOR_QUERY = "world/region{region@name}/sector[+NamedFilter,StringRegexMatches,_W$]/subsector{subsector@name}/technology{technology@name}/period/input[+NamedFilter,StringRegexMatches,water withdrawals$]/physical-demand{year@year}"
    WITHDRAWAL_SOURCE_QUERY = "world/region{region@name}/resource[+NamedFilter,StringRegexMatches,_water withdrawals]/subresource{subresource@name}//output/physical-output{year@year}"

    # f_maxsubresource = '/qfs/people/zhao924/gcam_water_storage/input/runoff_impacts_MIROC-ESM-CHEM_rcp6p0.csv'

    global gcam

    def __init__(self, period_n, supply_curve_lp, maxsubresource, save_outputs=None):

        self.period_n = period_n
        self.supply_curve_lp = supply_curve_lp
        self.save_outputs = save_outputs

        # get year
        self.year = 2020 + 5 * (self.period_n - 5)

        # self.maxsubresource_climate = self.load_maxsubresource(InfoExchange.f_maxsubresource)
        self.maxsubresource_climate = maxsubresource.loc[maxsubresource['year'].isin([self.year])]

        # run GCAM at each time period
        print('Running default GCAM for period: ', self.year)
        gcam.run_to_period(self.period_n)

        # water cost
        self.supply_curve_cost = gcam.get_data(InfoExchange.SUPPLY_COST_QUERY,
                                               **{"region": ["*"], "grade": ["*"]})

        # water available
        self.supply_curve_avail = gcam.get_data(InfoExchange.SUPPLY_AVAIL_QUERY,
                                                **{"region": ["*"], "grade": ["*"]})

        # maxSubResource
        self.maxsubresource = gcam.get_data(InfoExchange.MAXSUBRESOURCE_QUERY,
                                            **{"region": ["*"], "year": ["=", self.year]})

        # price
        self.price = gcam.get_data(InfoExchange.PRICE_SOLVE_QUERY,
                                   **{"year": ["=", self.year]})

        # initialization
        self.wd_sector_update = None
        self.wd_source_update = None
        self.price_update = None
        self.wd_sector_format = None
        self.wd_source_format = None
        
        # proceed with updating GCAM
        self.update_gcam()
        
        # get updated GCAM values
        self.get_gcam_update()
        
        # reformat demand output
        self.demand_gcam = self.format_demand()

    def load_maxsubresource(self, fn, header_num=4):
        """
        Load and format xanthos runoff to required format for maxsubresource by GCAM.

        :return:                dataframe for maxsubresource
        """
        if not os.path.isfile(fn):
            raise IOError("Error: File does not exist:", fn)

        # for CSV files
        elif fn.endswith('.csv'):
            df = pd.read_csv(fn, skiprows=header_num)
            df = df.rename(columns={'renewresource': 'resource',
                                    'sub.renewable.resource': 'subresource',
                                    'year.fillout': 'year'})
            df = df.loc[df['year'].isin([self.year])]

        else:
            raise RuntimeError("File {} has unrecognized extension".format(fn))

        return df

    def update_gcam(self):
        """
        Pass info to GCAM using gcamwrapper.
        """
        # update cost
        self.supply_curve_cost.set_index(['resource', 'grade'], inplace=True)
        self.supply_curve_cost.update(self.supply_curve_lp.set_index(['resource', 'grade']))
        self.supply_curve_cost.reset_index(inplace=True)
        self.supply_curve_cost = self.supply_curve_cost[['region', 'resource', 'subresource', 'grade', 'extractioncost']]
        gcam.set_data(self.supply_curve_cost, InfoExchange.SUPPLY_COST_QUERY,
                      **{"region": ["+", "="], "grade": ["+", "="]})

        # update available
        self.supply_curve_avail.set_index(['resource', 'grade'], inplace=True)
        self.supply_curve_avail.update(self.supply_curve_lp.set_index(['resource', 'grade']))
        self.supply_curve_avail.reset_index(inplace=True)
        self.supply_curve_avail = self.supply_curve_avail[['region', 'resource', 'subresource', 'grade', 'available']]
        gcam.set_data(self.supply_curve_avail, InfoExchange.SUPPLY_AVAIL_QUERY,
                      **{"region": ["+", "="], "grade": ["+", "="]})

        # update maxSubResource
        self.maxsubresource.set_index(['resource', 'year'], inplace=True)
        self.maxsubresource.update(self.maxsubresource_climate.set_index(['resource', 'year']))
        self.maxsubresource.reset_index(inplace=True)
        self.maxsubresource = self.maxsubresource[['region', 'resource', 'subresource', 'year', 'maxSubResource']]
        gcam.set_data(self.maxsubresource, InfoExchange.MAXSUBRESOURCE_QUERY,
                      **{"region": ["+", "="], "year": ["+", "="]})

        print('Running updated GCAM for period: ', self.year)
        gcam.run_to_period(self.period_n)

    def get_gcam_update(self):
        """
        Get solved demand after updating the GCAM info.
        :return:
        """

        # get sectoral water withdrawal
        self.wd_sector_update = gcam.get_data(InfoExchange.WITHDRAWAL_SECTOR_QUERY,
                                              **{"region": ["*"], "subsector": ["*"],
                                                 "technology": ["*"], "year": ["=", self.year]})
        self.wd_sector_update.to_csv(os.path.join(self.save_outputs, 'gcam_demand_sector_' + str(self.year) + '.csv'),
                                     encoding='utf-8', index=False)

        # get groundwater vs. runoff withdrawal
        self.wd_source_update = gcam.get_data(InfoExchange.WITHDRAWAL_SOURCE_QUERY,
                                              **{"region": ["*"], "subresource": ["*"],
                                                 "year": ["=", self.year]})
        self.wd_source_update.to_csv(os.path.join(self.save_outputs, 'gcam_demand_source_' + str(self.year) + '.csv'),
                                     encoding='utf-8', index=False)

        # get price
        self.price_update = gcam.get_data(InfoExchange.PRICE_SOLVE_QUERY,
                                          **{"year": ["=", self.year]})

        if self.save_outputs is not None:
            out = pd.merge(self.price_update, self.price, how='left',
                           on=('market', 'year'), suffixes=('.original', '.new'))
            out.to_csv(os.path.join(self.save_outputs, 'price_' + str(self.year) + '.csv'), encoding='utf-8',
                       index=False)

    def format_demand(self):
        """
        Merge sectoral water demand and source (groundwater vs runoff) demand.
        
        :return:                dataframe
        """
        # get solved sectoral water withdrawal from GCAM
        # note: irrigation withdrawal is by basin; other sector withdrawal is by region        
        self.wd_sector_format = self.wd_sector_update.loc[self.wd_sector_update['technology'] == 'water withdrawals']

        # get solved groundwater vs. runoff withdrawal from GCAM
        # use this withdrawal to get only surfacewater withdrawals
        self.wd_source_format = self.wd_source_update.loc[self.wd_source_update['subresource'] == 'runoff']
        self.wd_source_format[['subsector', 'technology']] = self.wd_source_format.resource.str.split(pat='_', expand=True).astype(str)
        self.wd_source_format = self.wd_source_format[['subsector', 'year', 'physical-output']]
        self.wd_source_format = self.wd_source_format.rename(columns={'physical-output': 'runoff'})

        # create demand from GCAM including sectoral demand and runoff withdrawals
        df = self.wd_sector_format.merge(self.wd_source_format, on=['subsector', 'year'])

        return df


def update_storage_capacity(basin_id, capacity_yield_curve, wd_source):
    """
    Calculate potential storage capacity expansion based on GCAM sovled water demand.

    :param basin_id:                        integer for basin id
    :param capacity_yield_curve:            dataframe for capacity yield curve
    :param wd_source:                       dataframe for GCAM solved sectoral water demand and source water demand

    :return:                                dataframe with basin id and storage capacity
    """

    # construct interpolatable yield-capacity curve
    yield_arr = capacity_yield_curve.loc[capacity_yield_curve['basin_id'] == basin_id]['yield'].values
    capacity_arr = capacity_yield_curve.loc[capacity_yield_curve['basin_id'] == basin_id]['capacity'].values
    yield_capacity_basin_interpld = interp1d(yield_arr, capacity_arr)

    # get gcam basin name
    gcam_basin_name = basin_mapping.loc[basin_mapping['basin_id'] == basin_id, ['gcam_basin_name']].iloc[0]['gcam_basin_name']

    # missing basins in GCAM output
    missing_basins = ['Guam', 'Micronesia', 'Antarctica']

    # get total surface water demand
    if gcam_basin_name not in missing_basins:
        # total surface water demand within the basin
        demand_basin = wd_source.loc[wd_source['subsector'] == gcam_basin_name]['runoff'].values[0]
        print('Min and Max yield for basin ', basin_id, 'is: ', yield_arr[0], ' and ', yield_arr[-1],
              '. GCAM solved demand from surface water is:', demand_basin)

        # determine storage capacity based on GCAM basin-level demand 
        if (demand_basin >= yield_arr[0]) & (demand_basin <= yield_arr[-1]):
            capacity_basin = yield_capacity_basin_interpld(demand_basin)
        elif demand_basin > yield_arr[-1]:
            capacity_basin = capacity_arr[-1]
        else:
            capacity_basin = 0
    else:
        demand_basin = 0
        capacity_basin = 0

    # construct dataframe
    df = pd.DataFrame({'basin_id': [basin_id],
                       'capacity_gcam': [capacity_basin],
                       'demand_gcam': [demand_basin],
                       'yield_min': [yield_arr[0]],
                       'yield_max': [yield_arr[-1]]})

    return df


if __name__ == '__main__':

    # initialize mpi4py parallelization
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()  # rank of processes from 0 to size-1
    size = comm.Get_size()  # total number of processes

    NUM_BASINS = 236
    block_size = int(NUM_BASINS / (size - 1))
    pd.set_option('display.max_columns', 6)

    feedback = True

    if rank != 0:
        start = (rank - 1) * block_size + 1
        if rank == size - 1:
            end = NUM_BASINS
        else:
            end = start + block_size

        # initialize variables
        demand_gcam = None
        capacity_gcam = None
        capacity_gcam_pre = None

        
    # initialize the output paths
    today = str(datetime.now().strftime('%Y%m%d_%Hh%Mm'))
    lp_path = os.path.join('your/path/to/outputs_glory/lp_outputs', today)
    gcam_path = os.path.join('your/path/to/outputs_glory/gcam_outputs', today)
    fig_path = os.path.join('your/path/to/outputs_glory/figure_outputs', today)

    # Run the model up to period 5, the year 2020, and loop until period 11 (2050)
    if rank == 0:

        # create output path if not existed
        if not os.path.exists(lp_path):
            os.makedirs(lp_path)
        if not os.path.exists(gcam_path):
            os.makedirs(gcam_path)
            
        for y in range(2020, 2055, 5):
            fig_path_period = os.path.join(fig_path, str(y))
            if not os.path.exists(fig_path_period):
                os.makedirs(fig_path_period)

        # Load basin mapping information
        basin_mapping = DataLoader.load_basin_mapping(DataLoader.f_basin_mapping)

        # Create a GCAM instance by giving it a configuration file and the appropriate working directory
        gcam = gcamwrapper.Gcam("configuration_ref.xml", "your/path/to/gcam-core/exe")

    for period_n in range(5, 12, 1):
        period = 2020 + 5 * (period_n - 5)

        supply_curve_block = np.empty((0, 4), dtype=float)
        capacity_yield_block = np.empty((3, 0), dtype=float)
        demand_profile_block = np.empty((0, 3), dtype=float)
        maxsubresource_block = np.empty((0, 5), dtype=float)
        cost_block = np.empty((0, 2), dtype=float)

        # run parallel on difference processes
        if rank != 0:
            print('This is rank:', rank, 'processing basins from: ', start, ' to ', end - 1)
            for id in range(start, end, 1):

                # create SupplyCurve object
                if feedback:
                    s = SupplyCurve(basin_id=id, period=period, demand_gcam=demand_gcam, capacity_gcam=capacity_gcam)
                else:
                    s = SupplyCurve(basin_id=id, period=period, demand_gcam=None, capacity_gcam=None)

                # get supply curves
                supply_curve_append = s.supply_curve # 2D array
                supply_curve_block = np.concatenate((supply_curve_block, supply_curve_append), axis=0)

                # get capacity-yield curves
                capacity_yield_append = np.array([s.expansion_seq, s.capacity_yield_interpld(s.expansion_seq), np.repeat(id, len(s.expansion_seq))])
                capacity_yield_block = np.concatenate((capacity_yield_block, capacity_yield_append), axis=1)
                
                # get demand profile
                demand_profile_append = np.array(list(s.d.demand_profile.items()))
                demand_profile_append = np.hstack((demand_profile_append, np.tile([[id]], 12).T))
                demand_profile_block = np.concatenate((demand_profile_block, demand_profile_append), axis=0)

                # get maxsubresource
                maxsubresource_append = s.maxsubresource
                if not maxsubresource_append[0][0] != maxsubresource_append[0][0]:
                    maxsubresource_block = np.concatenate((maxsubresource_block, maxsubresource_append), axis=0)
                    
                # get unit cost 1975$/m3
                cost_append = s.cost_per_m3
                cost_block = np.concatenate((cost_block, cost_append), axis=0)

            # record the storage capacity at (t-1)
            print('This is rank ', rank, ', saving GCAM solved capacity for time period: ', period)
            capacity_gcam_pre = capacity_gcam
            # if (period > 2020) and (capacity_gcam_pre is not None):
            #     print('This is rank ', rank, ', GCAM solved capacity for period ', period, ' is:',
            #           capacity_gcam_pre.loc[(capacity_gcam_pre['capacity_gcam'] > 0) &
            #                                 (capacity_gcam_pre['basin_id'] >= start) &
            #                                 (capacity_gcam_pre['basin_id'] < end)].head(10))
            # else:
            #     print('This is period 2020 or feedback is set to False, no GCAM solved capacity being saved.')

        # barrier to wait until all processes are completed
        comm.barrier()

        # gather results from all processes (output data type is list)
        supply_curve_list = comm.gather(supply_curve_block, root=0)
        capacity_yield_list = comm.gather(capacity_yield_block.T, root=0)
        demand_profile_list = comm.gather(demand_profile_block, root=0)
        maxsubresource_list = comm.gather(maxsubresource_block, root=0)
        cost_list = comm.gather(cost_block, root=0)

        # perform GCAM run on the rank=0 process
        if rank == 0:

            # format supply curve as dataframe
            supply_curve_copy = copy.deepcopy(supply_curve_list)
            supply_curve_df = pd.DataFrame(np.vstack(supply_curve_copy),
                                           columns=['resource', 'grade', 'available', 'extractioncost'])
            supply_curve_df['available'] = supply_curve_df['available'].astype('float64', errors = 'raise')
            supply_curve_df['extractioncost'] = supply_curve_df['extractioncost'].astype('float64', errors = 'raise')
            supply_curve_df.to_csv(
                os.path.join(lp_path, 'LP_supply_curve_' + str(period) + '.csv'),
                encoding='utf-8', index=False)
            print(supply_curve_df.head(40))

            # format capacity-yield curve as dataframe
            capacity_yield_copy = copy.deepcopy(capacity_yield_list)
            capacity_yield_df = pd.DataFrame(np.vstack(capacity_yield_copy),
                                             columns=['capacity', 'yield', 'basin_id'])
            capacity_yield_df.to_csv(os.path.join(lp_path, 'LP_capacity_yield_curve_' + str(period) + '.csv'),
                                     encoding='utf-8', index=False)
            print(capacity_yield_df.tail(5))
            
            # format demand profile as dataframe
            demand_profile_copy = copy.deepcopy(demand_profile_list)
            demand_profile_df = pd.DataFrame(np.vstack(demand_profile_copy),
                                             columns=['month', 'fraction', 'basin_id'])
            demand_profile_df.to_csv(os.path.join(lp_path, 'LP_demand_profile_' + str(period) + '.csv'),
                                     encoding='utf-8', index=False)
            print(demand_profile_df.head(5))

            # format maxsubresource
            maxsubresource_copy = copy.deepcopy(maxsubresource_list)
            maxsubresource_df = pd.DataFrame(np.vstack(maxsubresource_copy),
                                             columns=['region', 'resource', 'subresource', 'year', 'maxSubResource'])
            maxsubresource_df['year'] = maxsubresource_df['year'].astype('int64', errors = 'raise')
            maxsubresource_df['maxSubResource'] = maxsubresource_df['maxSubResource'].astype('float64', errors = 'raise')
            maxsubresource_df.to_csv(os.path.join(lp_path, 'LP_maxsubresource_' + str(period) + '.csv'),
                                     encoding='utf-8', index=False)
            print(maxsubresource_df.head(5))
            
            # format unit cost for reservoir construction
            cost_copy = copy.deepcopy(cost_list)
            cost_df = pd.DataFrame(np.vstack(cost_copy),
                                   columns=['basin_id', 'unit_cost'])
            cost_df['unit_cost'] = cost_df['unit_cost'].astype('float64', errors='raise')
            cost_df.to_csv(os.path.join(lp_path, 'LP_unit_cost_' + str(period) + '.csv'),
                           encoding='utf-8', index=False)

            # start information exchange with GCAM
            print('Starting gcamwrapper update ...........................')
            w = InfoExchange(period_n=period_n,
                             supply_curve_lp=supply_curve_df,
                             maxsubresource=maxsubresource_df,
                             save_outputs=gcam_path)

            # format demand from GCAM including sectoral demand and runoff withdrawals
            demand_gcam = w.demand_gcam

            # get formatted source demand
            wd_source = w.wd_source_format

            # initialize capacity dataframe 
            capacity_gcam = pd.DataFrame(columns=['basin_id', 'capacity_gcam', 'demand_gcam', 'yield_min', 'yield_max'])

            # get the storage capacity value based on the solved withdrawal from GCAM
            # to determine if there will be expansion
            for id in range(1, 236, 1):

                # get capacity from surface water demand
                capacity_append = update_storage_capacity(basin_id=id,
                                                          capacity_yield_curve=capacity_yield_df,
                                                          wd_source=wd_source)
                capacity_gcam = pd.concat([capacity_gcam, capacity_append])

            # save capacity from each time period
            capacity_gcam.to_csv(os.path.join(lp_path, 'LP_gcam_capacity_' + str(period) + '.csv'),
                                 encoding='utf-8', index=False)

        if feedback:
            # broadcast demand and capacity to other processes
            demand_gcam = comm.bcast(demand_gcam, root=0)
            capacity_gcam = comm.bcast(capacity_gcam, root=0)

    # save GCAM database
    if rank == 0:
        print('Writing GCAM output xml database ...........................')
        gcam.print_xmldb()


