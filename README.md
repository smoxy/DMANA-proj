https://www.r-econometrics.com/reproduction/wooldridge/wooldridge14/

http://www.ozanbakis.com/R_slides/R_panel.pdf

https://www.youtube.com/watch?v=ZEftL_fnNGQ

## Finite Mixture
- aggiungere gli altri modelli (tclust (trimming), teigen, pgmmEM)
- selezione variabili (kmeans sparse, clustvarsel, vscc, SelvarClustLasso, ELASTIC(?))
- creare funzione con dati e lista di modelli su cui iterare creando un df con righe nomi dei modelli e colonne randIndex, adjRandIndex, BIC, ICL
- conclusioni sulla base di: 1) fit dei dati con osservazioni; 2) criterio di selezione BIC o ICL
- risultato <- table(mario$Y, pippo$map)
- rownames(risultato) <- c('71','73','74','70','71','72','73','74','75','76','74','76','78','79')
- t(risultato)
