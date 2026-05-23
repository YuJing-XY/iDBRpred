/*
 * defs.h
 *
 *  Created on: Oct 31, 2015
 *      Author: nmalhis
 */

#ifndef DEFS_H_
#define DEFS_H_

#define F_SIZE	8
#define _AV 	0
#define _MX		1
#define _RBF	2
#define _Sigmoid	3
#define _SPDB	0
#define _URDB	1

#define PATH_DATA    "inData/"
#define PATH_TMP     "tmp/"
#define PATH_SVM     "svmModels/"
#define PATH_MC      "MC1/"


#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

// MCW
#define DistributionSIZE 10001


#endif /* DEFS_H_ */
