exports = module.exports = {};
const Matrix = require('ml-matrix').Matrix;
const LM = require('my-ml-curve-fitting');
const stat = require('statistics-distribution');

exports.locate = function (obj) {  
	/* obj you pass in should have a key 'data' and if you want weighting then also the key 'weights'; data sample: [ [1, 3, 4.1], [2, 4, 5], [3, 4, 5.1], [2, -2, 1]]
	data of format (double x, double y, double distance) , weights nothing or array, for 1/x weights for data sample above: [1/4.1,1/5,1/5.1,1/1]
	-> Returns x: x-Value of best estimated point, y: y-Value of best fit point, 
	x_confidence_interval: 95 % confidence Intervall for x value, y_confidence_interval: 95 % confidence Intervall for x value, chi_squared: for Solution parameters
	log: function for logging results
	*/
	var data;
	var weights;
	if (obj.data != undefined){
		data = obj.data;
	}
	else {
		if (Array.isArray(obj[0])) {
			data = obj[0];
			if (Array.isArray(obj[1])) weights = obj[1];
		}
		else {
			throw new Error("You did not pass in Data! The Object you pass in needs a Key named 'data' that contains an Array with your Data!");
		}
	}
	if (obj.weights != undefined) {
		weights = obj.weights;
	}


	//Distance function. p is the guessed/initial point.
	var euclidean = function(t,p,c){
	    var rows = t.rows;
	    var result = new Matrix(t.rows, 1);
	    for(var i=0;i<rows;i++){
	       result[i][0] = Math.sqrt(Math.pow(t[i][0]-p[0][0],2)+Math.pow(t[i][1]-p[1][0],2));
	    }
	    return result;
	};

	var nbPoints = data.length;
	if (nbPoints < 3) throw new Error(nbPoints + " are to few Data Points to Trilaterte!");
	var t = new Matrix(nbPoints,2);          // number of independent variables -> 2 in case of 2-Dimensional Trilateration

	var y_data = new Matrix(nbPoints,1); 	// Only y_data
	var weight= new Matrix(nbPoints,1);
	if (weights != undefined) {
		if (weights.length != nbPoints) throw new Error("Weights need to be same length as the data passed!");
		for(var i=0;i<weights.length;i++){
			weight[i][0] = weights[i];
		}
	}
	else {
		for(var i=0;i<nbPoints;i++){
	  	weight[i][0]=1;				// 1 for non weigthed and '1/data[i][2]' for weigthed (inversly proportional to distance/radius)
		}
	}

	var sumX = 0;
	var sumY = 0;
	for(var i=0;i<nbPoints;i++){
	    t[i][0] = data[i][0];
	    sumX += data[i][0];
	    t[i][1] = data[i][1];
	    sumY += data[i][1];
	    y_data[i][0]=data[i][2];
	}
	
	/*
	opts   = vector of algorithmic parameters
	parameter    				defaults   meaning
	 opts(1)  =  prnt            3        >1 intermediate results; >2 plots
	 opts(2)  =  MaxIter      10*Npar     maximum number of iterations
	 opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
	 opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
	 opts(5)  =  epsilon_3       1e-3     convergence tolerance for Chi-square
	 opts(6)  =  epsilon_4       1e-2     determines acceptance of a L-M step
	 opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
	 opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
	 opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
	 opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda updat
	*/
	var opts = [ 2, 500, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];
	var consts = [];
	var p_init = new Matrix([[(sumX/nbPoints)],[(sumY/nbPoints)]]);	//Initial guess of point, start for Levenberg–Marquardt algorithm; to get an OK initial guess we just use the average coordinates of the Beacons we triliterate around
	var p_min = new Matrix([[-1000000],[-1000000]]);	//Bounds for the Parameters to be found
	var p_max = new Matrix([[1000000],[1000000]]);
	var fit = LM.optimize(euclidean,p_init,t,y_data,weight,0.001,p_min,p_max,consts,opts);	//OPTIMIZE CALL! Algorithm start
	var fit_x2 = fit.X2;						//Chi Squared
	var fit_p = fit.p;							//Parameter of the best Fit -> Coordinates of best solution
	//var covariance_matrix = fit.covar;		//Covariance Matrix of solution
	var standardErrorX = fit.sigma_p[0];		//Asymptotic standard error of parameter x
	var standardErrorY = fit.sigma_p[1];		//Asymptotic standard error of parameter y
	
	var t_Value = stat.tdistr(nbPoints-2,.025);		//T-Value for n-2 Degrees of freedom (we generate 2 Parameters) at with 5% (2*2,5% - Both sides) Margin of Error		 
	var x_Range = [(Number(fit_p[0]) - (t_Value*standardErrorX)),(Number(fit_p[0]) + (t_Value*standardErrorX))]			// +/- T-Value(for n-2 degrees of freedom (95% accuracy))
	var y_Range = [(Number(fit_p[1]) - (t_Value*standardErrorY)),(Number(fit_p[1]) + (t_Value*standardErrorY))]

	var result = {x:fit_p[0], y:fit_p[1],x_confidence_interval:x_Range,y_confidence_interval:y_Range,chi_squared:fit_x2, log:function(){
		console.log('Trilateration Result:\nEstimated Point: \n x= ' + result.x + ' , y= ' + result.y + '\n chi_squared= ' + result.chi_squared);
		console.log(' x confidence Interval (95%): [' + result.x_confidence_interval[0] + ',' + result.x_confidence_interval[1] + ']');
		console.log(' y confidence Interval (95%): [' + result.y_confidence_interval[0] + ',' + result.y_confidence_interval[1] + ']');
	}};
	return result;
}

exports.locateEllipse = function (obj) {  
	/* obj you pass in should have a key 'data' and if you want weighting then also the key 'weights'
	data[i][0] = beacon Position x value	
	data[i][1] = beacon Position y value	
	data[i][2] = rotation anti clockwise in degrees
	data[i][3] = beam_width
	data[i][4] = length of ellipse
	weights nothing or array of absolute number of weights
	-> Returns x: x-Value of best estimated point, y: y-Value of best fit point, 
	x_confidence_interval: 95 % confidence Intervall for x value, y_confidence_interval: 95 % confidence Intervall for x value, chi_squared: for Solution parameters
	log: function for logging results
	*/
	var data;
	var weights;
	if (obj.data != undefined){
		data = obj.data;
	}
	else {
		if (Array.isArray(obj) && Array.isArray(obj[0])) {
			data = obj[0];
			if (Array.isArray(obj[1])) weights = obj[1];
		}
		else {
			throw new Error("You did not pass in Data! The Object you pass in needs a Key named 'data' that contains an Array with your Data!");
		}
	}
	if (obj.weights != undefined) {
		weights = obj.weights;
	}


	//Distance function. p is the guessed/initial point.
		var euclideanToEllipse = function(t,p,c){			// 	Function to fit to our Data
	    //t[i][0] = beaconPosX									see https://math.stackexchange.com/questions/2256297/minimum-and-maximum-distance-of-a-point-from-an-ellipse#answer-2256352
	    //t[i][1] = beaconPosY									and https://math.stackexchange.com/questions/426150/what-is-the-general-equation-of-the-ellipse-that-is-not-in-the-origin-and-rotate#answer-434482
	    //t[i][2] = rotation
	    //t[i][3] = beam_width
	    //t[i][4] = distance
	    //p[0][0] = pointX
	    //p[1][0] = pointY
	    var rows = t.rows;
	    var result = new Matrix(t.rows, 1);
		    for(var i=0;i<rows;i++){
		    	var b_w_h_r	= toRadians(t[i][3]/2); 	// half of Beam with in Radians
				var ax = t[i][4]/2;					//http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=4753864
				var ay = ax*Math.sqrt(
					Math.tan(b_w_h_r) * (
						Math.sin(b_w_h_r)/(
						Math.sqrt(2) - Math.cos(b_w_h_r)
						)
					)
				);
				var cosResult = Math.cos(toRadians(t[i][2]));		
				var sinResult = Math.sin(toRadians(t[i][2]));
				var pX = (((p[0][0]-t[i][0])*cosResult + (p[1][0]-t[i][1])*sinResult) - ax) / ax;
				var pY = ((p[1][0]-t[i][1])*cosResult - (p[0][0]-t[i][0])*sinResult) / ay;
				var magV = Math.sqrt(pX*pX + pY*pY);
				var pointOnCirlceX = pX / magV ;
				var pointOnCircleY = pY / magV ;
				var les = new Matrix (2,3); 	// Linear equatation system;
				les[0][0] = cosResult;								// Line 1: cosResult * x + sinResult * y = ...
				les[0][1] = sinResult;
				les[0][2] = (pointOnCirlceX*ax) + ax + (t[i][0]*cosResult) + (t[i][1]*sinResult);
				les[1][0] = -sinResult;								// Line 2 : -sinResult * x + cosResult * y = ...
				les[1][1] = cosResult;
				les[1][2] = (pointOnCircleY*ay) + (t[i][1]*cosResult) - (t[i][0] * sinResult);
				var factor = les[1][0]/les[0][0];	
				les.setRow(1,les.selection([1],[0,1,2]).apply(function (i,j) {
					this[i][j] -= factor * les[0][j] ;
				}));
				les[1][2] = les[1][2]/les[1][1];
				les[1][1] = 1;
				var yInOldCoordSystem = les[1][2];
				les[0][2] = (les[0][2] - (les[0][1] * les[1][2])) / les[0][0];
				les[0][1] = 0;
				les[0][0] = 1;
				var xInOldCoordSystem = les[0][2];
		       	result[i][0] = (Math.sqrt(Math.pow(xInOldCoordSystem-p[0][0],2)+Math.pow(yInOldCoordSystem-p[1][0],2)));
		    }
	    return result;
		};

	var nbPoints = data.length;
	if (nbPoints < 2) throw new Error(nbPoints + " are to few Data Points to Trilaterte!");
	var t = new Matrix(nbPoints,5);         

	var y_data = new Matrix(nbPoints,1);
	var weight= new Matrix(nbPoints,1);
	if (weights != undefined) {
		if (weights.length != nbPoints) throw new Error("Weights need to be same length as the data passed!");
		for(var i=0;i<weights.length;i++){
			weight[i][0] = weights[i];
		}
	}
	else {
		weight.fill(1);
	}

	var sumX = 0;
	var sumY = 0;
	for(var i=0;i<nbPoints;i++){
	    t[i][0] = data[i][0];
	    sumX += data[i][0];
	    t[i][1] = data[i][1];
	    sumY += data[i][1];
	    t[i][2] = data[i][2];
	    t[i][3] = data[i][3];
	    t[i][4] = data[i][4];
	    y_data[i][0]=0;
	}
	
	/*
	opts   = vector of algorithmic parameters
	parameter    				defaults   meaning
	 opts(1)  =  prnt            3        >1 intermediate results; >2 plots
	 opts(2)  =  MaxIter      10*Npar     maximum number of iterations
	 opts(3)  =  epsilon_1       1e-3     convergence tolerance for gradient
	 opts(4)  =  epsilon_2       1e-3     convergence tolerance for parameters
	 opts(5)  =  epsilon_3       1e-3     convergence tolerance for Chi-square
	 opts(6)  =  epsilon_4       1e-2     determines acceptance of a L-M step
	 opts(7)  =  lambda_0        1e-2     initial value of L-M paramter
	 opts(8)  =  lambda_UP_fac   11       factor for increasing lambda
	 opts(9)  =  lambda_DN_fac    9       factor for decreasing lambda
	 opts(10) =  Update_Type      1       1: Levenberg-Marquardt lambda updat
	*/
	var opts = [ 2, 500, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];
	var consts = [];
	var p_init = new Matrix([[(sumX/nbPoints)],[(sumY/nbPoints)]]);	//Initial guess of point, start for Levenberg–Marquardt algorithm; to get an OK initial guess we just use the average coordinates of the Beacons we triliterate around
	var p_min = new Matrix([[-1000000],[-1000000]]);	//Bounds for the Parameters to be found
	var p_max = new Matrix([[1000000],[1000000]]);
	var fit = LM.optimize(euclidean,p_init,t,y_data,weight,0.001,p_min,p_max,consts,opts);	//OPTIMIZE CALL! Algorithm start
	var fit_x2 = fit.X2;						//Chi Squared
	var fit_p = fit.p;							//Parameter of the best Fit -> Coordinates of best solution
	//var covariance_matrix = fit.covar;		//Covariance Matrix of solution
	var standardErrorX = fit.sigma_p[0];		//Asymptotic standard error of parameter x
	var standardErrorY = fit.sigma_p[1];		//Asymptotic standard error of parameter y
	
	var t_Value = stat.tdistr(nbPoints-2,.025);		//T-Value for n-2 Degrees of freedom (we generate 2 Parameters) at with 5% (2*2,5% - Both sides) Margin of Error		 
	var x_Range = [(Number(fit_p[0]) - (t_Value*standardErrorX)),(Number(fit_p[0]) + (t_Value*standardErrorX))]			// +/- T-Value(for n-2 degrees of freedom (95% accuracy))
	var y_Range = [(Number(fit_p[1]) - (t_Value*standardErrorY)),(Number(fit_p[1]) + (t_Value*standardErrorY))]

	var result = {x:fit_p[0], y:fit_p[1],x_confidence_interval:x_Range,y_confidence_interval:y_Range,chi_squared:fit_x2, log:function(){
		console.log('Trilateration Result:\nEstimated Point: \n x= ' + result.x + ' , y= ' + result.y + '\n chi_squared= ' + result.chi_squared);
		console.log(' x confidence Interval (95%): [' + result.x_confidence_interval[0] + ',' + result.x_confidence_interval[1] + ']');
		console.log(' y confidence Interval (95%): [' + result.y_confidence_interval[0] + ',' + result.y_confidence_interval[1] + ']');
	}};
	return result;
}


exports.learn = function (arr) {
	/*	Input is an array: For each of the Beacons that were used in the Trilateration 1 Object.
	This Object has a 'data' Property which is an array that contains all data points. Each data point has at index [0] the rssi data.
		At index [1] is the REAL distance between the beacon and the to-locate-object
 	This Object should also have a 'distanceFormula' Property which is an object that has the Propertys 'expDivisor' and 'divisor' ---> (10^(Abs(RSSI)/expDivisor))/divisor
 	This Object can also have a 'id' Property which will just be given trough, So you can identify later assign the results to a specific Beacon.
 	This Object can also have a 'weights' property which is an array that contains the according weights for the 'data' points. So it has to be same length as this array.

 	This function will return an array of the same length as the input array: For each Beacon it will calculate the optimal RSSI->Distance formula Parameters (expDivisor & divisor).
 	Each Object in the array has the propertys 'id', 'expDivisor', 'divisor'
	*/

	//Distance Formula. p is the guessed/initial point.
	var distanceFormula = function(t,p,consts){
	    var rows = t.rows;
	    var result = new Matrix(t.rows, 1);
	    for(var i=0;i<rows;i++){
	      result[i][0] = Math.pow(10,(Math.abs(t[i][0])/p[0][0]))/p[1][0];
	    }
	    return result;
	};

	var arr_length = arr.length;
	var results = [];

	for (var w=0; w< arr_length;w++) {
		if (arr[w].hasOwnProperty("data")) {
			var data = arr[w].data;
		}
		else throw new Error("No 'data' property in input!");

		var data_length = data.length;
		var t = new Matrix(data_length,1);     
		var y_data = new Matrix(data_length,1)
		var weight= new Matrix(data_length,1);
		var weights;
		if (arr[w].weights != undefined) {
			weights = arr[w].weights;
			if (weights === undefined || weights.length != data_length) throw new Error("Weights need to be same length as the data passed!");
			for(var i=0;i<weights.length;i++){
				weight[i][0] = weights[i];
			}
		}
		else {
			for(var i=0;i<data_length;i++){
		  		weight[i][0]=1;				// 1 for non weigthed and '1/data[i][2]' for weigthed (inversly proportional to distance/radius)
			}
		}

		for(var i=0;i<data_length;i++){
		    t[i][0] = data[i][0];
		    y_data[i][0]=data[i][1];
		}
		var init = [];
		if (arr[w].hasOwnProperty("distanceFormula")) {
			init[0]=arr[w].distanceFormula.expDivisor;
			init[1]=arr[w].distanceFormula.divisor;
		}
		else {
			init[0]=20;
			init[1]=1000;
		}


		var opts = [ 2, 500, 1e-3, 1e-3, 1e-3, 1e-2, 1e-2, 11, 9, 1 ];
		var consts = [];
		var p_init = new Matrix([[init[0]],[init[1]]]);			//Initial guess of point, start for Levenberg–Marquardt algorithm
		var p_min = new Matrix([[0],[0]]);	//Bounds for the Parameters to be found
		var p_max = new Matrix([[100],[1000000]]);
		var fit = LM.optimize(distanceFormula,p_init,t,y_data,weight,0.001,p_min,p_max,consts,opts);	//OPTIMIZE CALL! Algorithm start
		//var fit_x2 = fit.X2;						//Chi Squared
		var fit_p = fit.p;							//Parameter of the best Fit 
		//var covariance_matrix = fit.covar;			//Covariance Matrix of solution
		//var standardErrorX = fit.sigma_p[0];		//Asymptotic standard error of parameter 
		//var standardErrorY = fit.sigma_p[1];		//Asymptotic standard error of parameter
		var id;
		if (arr[w].hasOwnProperty("id")) {
			id = arr[w].id;
		} 
		results[w]={id:id,expDivisor:fit_p[0],divisor:fit_p[1]};
	}	//End outer loop (number of Beacons to be adjusted)
	//console.log(results);

	return results;
}