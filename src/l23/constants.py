"""Constants shared across steps: Yoo25 L2/3 cells and ages, the NR vs DR test design, colours."""

# Yoo25 and Cheng22 L2/3 cells of the combined snRNA object (superdupermegaRNA_hasraw.h5ad)
YOO25_STUDY = '2023 Multiome'
CHENG22_AGE = 'P28'        # Cheng22 NR (Study missing); P28_dr / P28_dl are DR / dark-light
L23_SUBCLASS = 'L2/3'
# obs Age -> (postnatal day, condition); NR at P6-P21, DR only at P12-P21
AGES = {'P6': (6, 'NR'), 'P8': (8, 'NR'), 'P10': (10, 'NR'), 'P12': (12, 'NR'), 'P14': (14, 'NR'),
        'P17': (17, 'NR'), 'P21': (21, 'NR'),
        'P12DR': (12, 'DR'), 'P14DR': (14, 'DR'), 'P17DR': (17, 'DR'), 'P21DR': (21, 'DR')}

# NR vs DR over age (l23.timecourse): Ages with both NR and DR, and the centering day = mean DR day
TEST_DAYS = [12, 14, 17, 21]
D0 = 16
FDR = 0.05
STARS = [(0.001, '***'), (0.01, '**'), (0.05, '*')]

COND_COLORS = {'NR': 'C1', 'DR': 'black'}
ARCH_COLORS = {'A': 'C0', 'B': 'C1', 'C': 'C2'}
