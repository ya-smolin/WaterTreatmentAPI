The current project is developed by the Ph.D. student in water treatment - Yaroslav Smolin. I study mathematical modeling and optimization of water treatment processes that allow to define the best purification settings and to reduce expenses on plant building. I concentrate on water treatment methods that can efficiently and economically remove toxic synthetic organic substances (biofiltration and adsorption).
The main idea is to share code with a community and store developed projects in a cloud.

#isothermAPI
It's a program with UI that makes work with adsorption isotherms very easy. It provides tools for finding the best isotherm type for the current process using different statistical characteristics. Program calculate the confidence interval, R-square, SSE. You can set aqua solubility parameter or set autofitting.

used: MATLAB + fitting and optimization toolbox

Support isoterms:
* BET
* Freundlich
* Freundlich-linear
* Langmure
* Freudlich-Langmure
* Langmure-Linear
* Linear
* Polanyi
* Polanyi-linear
* Toth
* Temkin
* Redlich-Peterson
* You can easily add your own isoterms, just modify Model.m file.

#Bioadsoption model
Implemented improved bioadsoption model using developments of Institute of Colloid and Water Chemistry.

#Biodegradetion parameters fitting
Implemented ideas from paper "4-nitrophenol biodegradation in a sequencing batch reactor:
kinetic study and effect of filling time" by M. Concetta Tomeia, M. Cristina Annesinib, S. Bussolettib
