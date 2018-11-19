import sys 
sys.path.insert(0,'../Stage_0')
from dict2obj      import obj 
from cost_class    import costs
#=======================================================================
# Set and run HYDRA
#=======================================================================

class _set_inputs:

   def first_init(self):

#====================================================================
# engine model initializations
#====================================================================

        self.initEngine()

#====================================================================
# cost analysis: initializations
#====================================================================

        ad                = self.all_dict 
        ops               = obj(ad['operations'])
        init_cost         = ad['purchase']
        beta_factors      = init_cost['Beta_acq_factors']
        self.costs        = costs(ops, init_cost, beta_factors, self.mission.cost_duration, \
        					self.rotor.nrotors)
        self.wt_redund    = ad['redund']

        return None

#=======================================================================
# END OF FILE
#=======================================================================