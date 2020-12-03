
# Motrpac colors abbreviations

assay_abbr = c(atacseq='ATAC',
               'atac-seq'='ATAC',
               rrbs='METHYL',
               rnaseq='TRNSCRPT',
               'rna-seq'='TRNSCRPT',
               'transcript-rna-seq'='TRNSCRPT',
               altsplice='SPLICE',
               metabolomics='METAB',
               'metab-meta'='METAB',
               'unnamed-metab'='U-METAB',
               'named-metab'='N-METAB',
               'prot-pr'='PROT',
               'prot-ac'='ACETYL',
               'prot-ph'='PHOSPHO')


save(assay_abbr, 
     file = "data/assay_abbr.RData" )

assay_order = c('METHYL',
                'ATAC',
                'SPLICE',
                'TRNSCRPT',
                'PROT',
                'PHOSPHO',
                'ACETYL',
                'N-METAB',
                'U-METAB')

save(assay_order, 
     file = "data/assay_order.RData" )

group_abbr = c('control'='SED',
               'control8w'='SED',
               '1w'='1W',
               '2w'='2W',
               '4w'='4W',
               '8w'='8W')

save(group_abbr, 
     file = "data/group_abbr.RData" )

group_cols = c('SED'='white',
               'control'='white',
               '0W' = 'white',
               '1W' = "#F7FCB9",
               '1w' = "#F7FCB9",
               '2W' = "#ADDD8E",
               '2w' = "#ADDD8E",
               '4W' = "#238443",
               '4w' = "#238443",
               '8W' = "#002612",
               '8w' = "#002612")

save(group_cols, 
     file = "data/group_cols.RData" )

sex_abbr = c('female'='F',
             'Female'='F',
             'male'='M',
             'Male'='M')

save(sex_abbr, 
     file = "data/sex_abbr.RData" )

sex_cols = c('M' = '#5555ff',
             'Male' = '#5555ff',
             'male' = '#5555ff',
             'F' = '#ff6eff',
             'Female' = '#ff6eff',
             'female' = '#ff6eff')

save(sex_cols, 
     file = "data/sex_cols.RData" )

tissue_abbr = c(adrenals='ADRNL',
                adrenal='ADRNL',
                aorta='VENACV',
                'vena-cava'='VENACV',
                'vena cava'='VENACV',
                'brown-adipose'='BAT',
                'brown adipose'='BAT',
                colon='COLON',
                cortex='CORTEX',
                plasma='PLASMA',
                'edta plasma'='PLASMA',
                gastrocnemius='SKM-GN',
                heart='HEART',
                hippocampus='HIPPOC',
                hypothalmus='HYPOTH',
                hypothalamus='HYPOTH',
                kidney='KIDNEY',
                liver='LIVER',
                lung='LUNG',
                ovaries='OVARY',
                'blood-rna'='BLOOD',
                'blood rna'='BLOOD',
                'paxgene rna'='BLOOD',
                'paxgene-rna'='BLOOD',
                'small-intestine'='SMLINT',
                'small intestine'='SMLINT',
                spleen='SPLEEN',
                testes='TESTES',
                spleen='SPLEEN',
                'vastus-lateralis'='SKM-VL',
                'vastus lateralis'='SKM-VL',
                'white-adipose'='WAT-SC',
                'white adipose'='WAT-SC')

save(tissue_abbr, 
     file = "data/tissue_abbr.RData" )

tissue_cols = c('paxgene rna'='#FF0000', # blood and paxgene rna are the same (RNA-seq)
  'paxgene-rna'='#FF0000',
  'paxgene_rna'='#FF0000',
  blood='#FF0000',
  'blood rna'='#FF0000',
  'blood-rna'='#FF0000',
  'blood_rna'='#FF0000',
  BLOOD='#FF0000',
  'edta plasma'="#FF6161",
  'edta-plasma'="#FF6161",
  'edta_plasma'="#FF6161",
  plasma="#FF6161",
  PLASMA="#FF6161",
  heart='#EB7602',
  HEART='#EB7602',
  aorta="#B22222",
  'vena cava'="#B22222",
  'vena-cava'="#B22222",
  'vena_cava'="#B22222",
  VENACV="#B22222",
  spleen="#0000FF",
  SPLEEN="#0000FF",
  gastrocnemius='#00CD00',
  'SKM-GN'='#00CD00',
  'vastus lateralis'='#077D07',
  'vastus-lateralis'='#077D07',
  'vastus_lateralis'='#077D07',
  'SKM-VL'='#077D07',
  'white adipose'='#1E90FF',
  'white-adipose'='#1E90FF',
  'white_adipose'='#1E90FF',
  'WAT-SC'='#1E90FF',
  'brown adipose'='#FFA500',
  'brown-adipose'='#FFA500',
  'brown_adipose'='#FFA500',
  'BAT'='#FFA500',
  liver="#E349E3",
  LIVER="#E349E3",
  lung="#14AE9E",
  LUNG="#14AE9E",
  kidney='#A020F0',
  KIDNEY='#A020F0',
  adrenal="#DDA0DD",
  adrenals="#DDA0DD",
  ADRNL="#DDA0DD",
  cortex='#F2CA00',
  CORTEX='#F2CA00',
  hypothalamus='#C1A100',
  HYPOTH='#C1A100',
  hypothalmus='#C1A100',
  hippocampus='#937B00',
  HIPPOC='#937B00',
  'small intestine'='#A18277',
  'small-intestine'='#A18277',
  'small_intestine'='#A18277',
  'SMLINT'='#A18277',
  colon='#7D4C3B',
  COLON='#7D4C3B',
  ovaries="#FFAEB9",
  OVARY="#FFAEB9",
  testes="#000080",
  TESTES="#000080")

save(tissue_cols, 
     file = "data/tissue_cols.RData" )

tissue_order = c('BLOOD',
                 'PLASMA',
                 'HEART',
                 'VENACV',
                 'SPLEEN',
                 'SKM-GN',
                 'SKM-VL',
                 'WAT-SC',
                 'BAT',
                 'LIVER',
                 'LUNG',
                 'KIDNEY',
                 'ADRNL',
                 'CORTEX',
                 'HYPOTH',
                 'HIPPOC',
                 'SMLINT',
                 'COLON',
                 'OVARY',
                 'TESTES')

save(tissue_order, 
     file = "data/tissue_order.RData" )

