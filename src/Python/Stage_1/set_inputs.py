import sys 
sys.path.insert(0,'../Stage_0')
from dict2obj       import obj 
from cost_class     import costs
from fuselage_class import geometry 

#=======================================================================
# Set and run HYDRA
#=======================================================================

class _set_inputs:

   def first_init(self):

#====================================================================
# cost analysis: initializations
#====================================================================

        ad                = self.all_dict 
        ops               = obj(ad['operations'])
        npax              = ad['aircraft']['pax_count']
        init_cost         = ad['purchase']
        if(bool(init_cost)):
            beta_factors      = init_cost['Beta_acq_factors']
            self.costs        = costs(ops, init_cost, beta_factors, self.mission, \
            					self.rotor.nrotors, npax)
        else:
            self.costs        = {}

#====================================================================
# assign geometry information: rotor layout
#====================================================================

        if('Geometry' in ad['empirical']):
            self.geometry         = geometry(ad['empirical']['Geometry'], self.rotor.nrotors, self.wing.nwings)
        else:
            print('setting default geometry for fuselage: may not be used in sizing?')
            self.geometry         = obj({'fuselage_length':1.0, 'fuselage_width':0.0, 'clearance':0.1})

#====================================================================
# tech factors: default initializations
#====================================================================
        
        emp             = self.emp_data
        keys            = ['rotor','wing','empennage','fuselage','landing_gear',       \
                           'fuel_system' ,'drive_system','flight_control','anti_icing', \
                           'powerplant'  ,'fuel','battery','emergency_sys']

        data            = {} 
        for k in keys:
            data[k]     = 1.0 
        defaults        = {'Weight_scaling': data}
        defaults        = obj(defaults)

#if tech factors not defined in input file, use default of unity
        if not(hasattr(emp, 'Tech_factors')):
            emp.Tech_factors = defaults 
            print('setting all tech factors to default')
        else:
# if weight scaling is not defined in tech factors, use default of unity
            if( not(hasattr(emp.Tech_factors, 'Weight_scaling'))):
                emp.Tech_factors.Weight_scaling = defaults.Weight_scaling

# here, some elements may/may not be found 
            else:
                for k in keys:
                    if(not(hasattr(emp.Tech_factors.Weight_scaling,k))):
                        setattr(emp.Tech_factors.Weight_scaling,k,defaults.Weight_scaling.__getattribute__(k))
                        print('setting tech factor for weights = 1, component is ',k)

#====================================================================
# redundancies: default initializations
#====================================================================

        rdict             = ad['redund']
        keylist           = ['wing_flap', 'tilt_actuator', 'wires', 'avionics']
        
        for k in keylist:
          if k not in rdict:
            rdict[k]      = 1.0

        self.wt_redund    = rdict 

        return None

#=======================================================================
# END OF FILE
#=======================================================================