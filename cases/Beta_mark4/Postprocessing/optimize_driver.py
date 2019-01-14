#====================================================================
# quantify sensitivity of design to various uncertainties in modeling
#====================================================================

import sys,os
sys.path.append('../../src/Python/Postprocessing/')
sys.path.append('../../src/Python/Stage_0/')
sys.path.append('../../src/Python/Stage_1/')
sys.path.append('../../src/Python/Stage_1/afdd/')
sys.path.append('../../src/Python/Stage_1/engines/')
sys.path.append('../../src/Python/Stage_3/')
sys.path.append(os.getcwd())

from optimize import optimize

optimize(os.getcwd())

