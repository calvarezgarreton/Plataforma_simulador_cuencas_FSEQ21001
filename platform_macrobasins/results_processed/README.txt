README
25 Mayo 2023

El archivo data_platform_mon_LTM.csv contiene los resultados de las simulaciones hidrológicas de macrocuencas seleccionadas de la bbdd CAMELS-CL para la herramienta "SImulador de cuencas" del proyecto FSEQ210001.
La modelación se hizo con el modelo mHM.
Autores: C. Alvarez-Garreton, J.P. Boisier


Las columnas del archivos contienen la sgte. información:

-------------------
Datos de la cuenca

gauge_id: código de la cuenca CAMELS-CL 
gauge_name: nombre de la estación fluviométrica de la salida de la cuenca.
area_km2: área en km2

-------------------------------------
Información del experimento

sim_short_name: nombre del experimento
sim_pr_name: fuente de datos de la precipitación
sim_pet_name: fuente de datos de la evapotranspiración potencial
sim_txn_name: fuente de datos de la temperatura
sim_lai_name: fuente de datos del índice de área foliar

-------------------------------------
Flujo medio mensual a escala de cuenca, calculado para los siguientes períodos:
1960-2020
1960-1990
1990-2020
2040-2070
2070-2100

Notar que algunos experimentos no tienen datos en ciertos períodos (e.g., simulación histórica no tiene datos en los períodos futuros 2040-2070 y 2070-2100)

date: mes                                  
pr_mm_mon_LTM_xxxx_yyyy: precipitación media mensual período xxxx a yyyy
et_mm_mon_LTM_xxxx_yyyy: evapotranspiración media mensual período xxxx a yyyy
tmin_C_mon_LTM_xxxx_yyyy: temperatura mínima mensual período xxxx a yyyy
tmax_C_mon_LTM_xxxx_yyyy: temperatura máxima mensual período xxxx a yyyy
tavg_C_mon_LTM_xxxx_yyyy: temperatura promedio mensual período xxxx a yyyy
pet_mm_mon_LTM_xxxx_yyyy: evapotranspiración potencial media mensual período xxxx a yyyy
qobs_m3s_mon_LTM_xxxx_yyyy: caudal observado medio mensual período xxxx a yyyy
qobs_mm_mon_LTM_xxxx_yyyy: escorrentía observada media mensual período xxxx a yyyy
qsim_m3s_mon_LTM_xxxx_yyyy: caudal simulado medio mensual período xxxx a yyyy
qsim_mm_mon_LTM_xxxx_yyyy: escorrentía simulada media mensual período xxxx a yyyy
