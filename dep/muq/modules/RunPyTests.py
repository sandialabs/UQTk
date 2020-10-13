import sys
sys.path.append('modules/Modeling/test')

import unittest

from PDEModPieceTests import *

suite = unittest.TestLoader().loadTestsFromTestCase(PDEModPieceTests)
unittest.TextTestRunner(verbosity=2).run(suite)
