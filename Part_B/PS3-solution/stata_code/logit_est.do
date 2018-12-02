eststo clear 
eststo: logit default debt_to_gdp, vce(robust)
esttab using lotable1.tex, replace label  title(Logit Regression Table\label{tab1}) booktabs alignment(D{.}{.}{-1}) width(0.8\hsize)
