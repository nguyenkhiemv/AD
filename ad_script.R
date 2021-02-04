#read data
ad_dat = read.csv("HBSpring2020_2020-06-18.T1.WMH.v2.v6.SNI.Repro.moregene.meso_DEIDENTIFIED.csv")

library( dplyr )

#try imputation with bmi and GMV_TIV_out 
data <- dplyr::select( ad_dat, c( 'PIDN', 'age_rounded', 'GMV_TIV_out', 'bmi' ) ) %>% data.frame( )
data <- na.omit( data ) #if NA is not removed, model2 and model3 can not run 

library( fcomplete )
data$GMV_TIV_out = as.numeric(as.character(data$GMV_TIV_out))
data$age_rounded

model <- fregression( GMV_TIV_out  ~ age_rounded | PIDN, data = data , 
                      d = 6, K = 3, lambda= c(0.1,0.2), cv.ratio=0.1, bins = 80 )

model2 <- fregression( GMV_TIV_out + bmi  ~ age_rounded | PIDN, data = data , d = 6, K = 3, 
                       lambda= seq( 0, 1, by = 0.05), cv.ratio=0.2, verbose = 0, bins = 80 )

model3 <- fregression( GMV_TIV_out ~ age_rounded + U1 | PIDN, data = data, model2$u, d = 6, K = 3,
                        lambda= seq( 0, 1, by = 0.05), cv.ratio=0.2, verbose = 0, bins = 80 )
