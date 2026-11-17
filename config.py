"""All constants for the Z boson dimuon project."""
DATA_URL = "https://opendata.cern.ch/record/5201/files/Dimuon_DoubleMu.csv"
RECORD_URL = "https://opendata.cern.ch/record/5201"
DATA_FILE = "data/Dimuon_DoubleMu.csv"
PDG_MZ = 91.1879
PDG_MZ_ERR = 0.0020
PDG_GAMMA_Z = 2.4955
PDG_MJPSI = 3.096900
PDG_MJPSI_ERR = 0.000006
ETA_MAX = 2.4
# Z measurement
Z_PT_MIN = 20.0
Z_RANGE = (70.0, 110.0)
Z_BIN_WIDTH = 0.5
Z_RANGE_VARIANTS = [(75.0, 105.0), (66.0, 116.0)]
Z_PT_VARIANTS = [25.0, 30.0]
NORM_GRID_POINTS = 4001
N_TOYS_Z = 200
# J/psi scale check
JPSI_RANGE = (2.9, 3.3)
JPSI_BIN_WIDTH = 0.005
CB_N = 3.0
# bump scan
SCAN_RANGE = (15.0, 60.0)
SCAN_BIN_WIDTH = 0.5
SCAN_MASS_FIRST = 16.0
SCAN_MASS_LAST = 59.0
SCAN_MASS_STEP = 0.5
N_TOYS_SCAN = 500
RANDOM_SEED = 42
GOF_ALPHA = 0.01

# Further constants fixed by the design text (not free choices)
DETECTOR_REGION_ETA = 1.2          # variant 6: both muons with |eta| below this
MIN_EXPECTED_PER_BIN = 5           # Pearson chi-square bins need at least this expected count
GOF_POINTS_PER_BIN = 21            # trapezoid points per bin for the model integral
M_DIFF_TOLERANCE = 0.01            # GeV, file M against recomputed M
RETRY_COUNT = 4                    # download tries
RETRY_WAIT_S = 20
# Nominal Z fit start values and limits
Z_MASS_START, Z_MASS_LIMITS = 91.0, (88.0, 94.0)
Z_SIGMA_START, Z_SIGMA_LIMITS = 1.5, (0.2, 5.0)
Z_LAM_START, Z_LAM_LIMITS = -0.02, (-0.5, 0.5)
Z_SIG_FRAC_START, Z_BKG_FRAC_START = 0.9, 0.1
Z_YIELD_MAX_FACTOR = 2.0
Z_WIDTH_FREE_LIMITS = (0.5, 6.0)
LINEAR_C_LIMITS = (-1.0, 1.0)
SECOND_C1_LIMITS, SECOND_C2_LIMITS = (-0.9, 0.9), (-0.4, 0.4)
# J/psi fit start values and limits
JPSI_MEAN_START, JPSI_MEAN_LIMITS = 3.097, (3.0, 3.2)
JPSI_SIGMA_START, JPSI_SIGMA_LIMITS = 0.03, (0.005, 0.1)
JPSI_BETA_START, JPSI_BETA_LIMITS = 1.5, (0.3, 5.0)
# Scan background polynomial
SCAN_X_CENTRE, SCAN_X_HALF = 37.5, 22.5
SCAN_GOF_NDOF = 86
TOY_STREAM_Z = 1
TOY_STREAM_SCAN = 2
# J/psi background and yield starts (not given in the instructions, see RUN_LOG.md)
JPSI_SIG_FRAC_START, JPSI_BKG_FRAC_START = 0.7, 0.3
JPSI_LAM_START, JPSI_LAM_LIMITS = -1.0, (-20.0, 20.0)
# A parameter closer than this fraction of its allowed range to a limit is reported as near a limit
NEAR_LIMIT_FRACTION = 0.03
