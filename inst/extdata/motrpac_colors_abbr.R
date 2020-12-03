
# Motrpac colors abbreviations

assay_abbr = c(atac='ATAC',
               atacseq='ATAC',
               'atac-seq'='ATAC',
               
               rrbs='METHYL',
               
               rnaseq='TRNSCRPT',
               'rna-seq'='TRNSCRPT',
               'transcript-rna-seq'='TRNSCRPT',
               
               altsplice='SPLICE',
               
               metabolomics='METAB',
               'metab'='METAB',
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
               '0w'='SED',
               '1w'='1W',
               '2w'='2W',
               '4w'='4W',
               '8w'='8W')

save(group_abbr, 
     file = "data/group_abbr.RData" )

group_cols = c('SED'='white',
               'control'='white',
               '0w' = 'white',
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
                "t60-adrenal"='ADRNL',
                
                aorta='VENACV',
                'vena-cava'='VENACV',
                'vena cava'='VENACV',
                "t65-aorta"='VENACV',
                
                'brown-adipose'='BAT',
                'brown adipose'='BAT',
                "t69-brown-adipose"='BAT',
                
                colon='COLON',
                "t61-colon"='COLON',
                
                cortex='CORTEX',
                "t53-cortex"='CORTEX',
                
                plasma='PLASMA',
                'edta plasma'='PLASMA',
                "t31-plasma"='PLASMA',
                
                gastrocnemius='SKM-GN',
                "t55-gastrocnemius"='SKM-GN', 
                
                heart='HEART',
                "t58-heart"='HEART',
                
                hippocampus='HIPPOC',
                "t52-hippocampus"='HIPPOC',
                
                hypothalmus='HYPOTH',
                hypothalamus='HYPOTH',
                "t54-hypothalamus"='HYPOTH',
                
                kidney='KIDNEY',
                "t59-kidney"='KIDNEY',
                
                liver='LIVER',
                "t68-liver"='LIVER',
                
                lung='LUNG',
                "t66-lung"='LUNG',
                
                ovaries='OVARY',
                "t64-ovaries"='OVARY',
                
                'blood-rna'='BLOOD',
                'blood rna'='BLOOD',
                'paxgene rna'='BLOOD',
                'paxgene-rna'='BLOOD',
                "t30-blood-rna"='BLOOD',
                
                'small-intestine'='SMLINT',
                'small intestine'='SMLINT',
                "t67-small-intestine"='SMLINT',
                
                spleen='SPLEEN',
                "t62-spleen"='SPLEEN',
                
                testes='TESTES',
                "t63-testes"='TESTES',
                
                spleen='SPLEEN',
                "t62-spleen"='SPLEEN',
                
                'vastus-lateralis'='SKM-VL',
                'vastus lateralis'='SKM-VL',
                "t56-vastus-lateralis"='SKM-VL',
                
                'white-adipose'='WAT-SC',
                'white adipose'='WAT-SC',
                "t70-white-adipose"='WAT-SC')

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
                "t30-blood-rna"='#FF0000',
                
                'edta plasma'="#FF6161",
                'edta-plasma'="#FF6161",
                'edta_plasma'="#FF6161",
                plasma="#FF6161",
                PLASMA="#FF6161",
                "t31-plasma"="#FF6161", 
                
                heart='#EB7602',
                HEART='#EB7602',
                "t58-heart"='#EB7602',
                
                aorta="#B22222",
                'vena cava'="#B22222",
                'vena-cava'="#B22222",
                'vena_cava'="#B22222",
                VENACV="#B22222",
                "t65-aorta"="#B22222",
                
                spleen="#0000FF",
                SPLEEN="#0000FF",
                "t62-spleen"="#0000FF",
                
                gastrocnemius='#00CD00',
                'SKM-GN'='#00CD00',
                "t55-gastrocnemius"='#00CD00',
                
                'vastus lateralis'='#077D07',
                'vastus-lateralis'='#077D07',
                'vastus_lateralis'='#077D07',
                'SKM-VL'='#077D07',"t56-vastus-lateralis"='#077D07',
                
                'white adipose'='#1E90FF',
                'white-adipose'='#1E90FF',
                'white_adipose'='#1E90FF',
                'WAT-SC'='#1E90FF',
                "t70-white-adipose"='#1E90FF',
                
                'brown adipose'='#FFA500',
                'brown-adipose'='#FFA500',
                'brown_adipose'='#FFA500',
                'BAT'='#FFA500',
                "t69-brown-adipose"='#FFA500',
                
                liver="#E349E3",
                LIVER="#E349E3",
                "t68-liver"="#E349E3",
                
                lung="#14AE9E",
                LUNG="#14AE9E",
                "t66-lung"="#14AE9E",
                
                kidney='#A020F0',
                KIDNEY='#A020F0',
                "t59-kidney"='#A020F0', 
                
                adrenal="#DDA0DD",
                adrenals="#DDA0DD",
                ADRNL="#DDA0DD",
                "t60-adrenal"="#DDA0DD",
                
                cortex='#F2CA00',
                CORTEX='#F2CA00',
                "t53-cortex"='#F2CA00', 
                
                hypothalamus='#C1A100',
                HYPOTH='#C1A100',
                hypothalmus='#C1A100',
                "t54-hypothalamus"='#C1A100',
                
                hippocampus='#937B00',
                HIPPOC='#937B00',
                "t52-hippocampus"='#937B00',
                
                'small intestine'='#A18277',
                'small-intestine'='#A18277',
                'small_intestine'='#A18277',
                'SMLINT'='#A18277',
                "t67-small-intestine"='#A18277',
                
                colon='#7D4C3B',
                COLON='#7D4C3B',
                "t61-colon"='#7D4C3B',
                
                ovaries="#FFAEB9",
                OVARY="#FFAEB9",
                "t64-ovaries"="#FFAEB9",
                
                testes="#000080",
                TESTES="#000080",
               "t63-testes"="#000080")

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

