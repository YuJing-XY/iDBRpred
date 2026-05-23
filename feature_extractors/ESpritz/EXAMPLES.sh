############################### Sequence only (very fast) ###########################

#execute the following for short xray disorder predictions which maxmixe Sw measure on training set
./espritz.pl example_fastas/ X 1
#execute the following for short xray disorder predictions at a 5% False Positive Rate (FPR)
./espritz.pl example_fastas/ X 0
#execute the following for long disprot type disorder predictions which maxmixe Sw measure on training set
./espritz.pl example_fastas/ D 1
#execute the following for long disprot type disorder predictions at a 5% FPR
./espritz.pl example_fastas/ D 0
#execute the following for NMR disorder predictions which maxmixe Sw measure on training set
./espritz.pl example_fastas/ N 1
#execute the following for NMR disorder predictions at a 5% FPR
./espritz.pl example_fastas/ N 0




############################### psi-blast (Slightly more accurate but slower) ###########################
#execute the following for short xray disorder predictions using psi-blast  which maxmixe Sw measure on training set
./espritz.pl example_fastas/ pX 1
#execute the following for short xray disorder predictions using psi-blast at a 5% False Positive Rate (FPR)
./espritz.pl example_fastas/ pX 0
#execute the following for long disprot type disorder predictions using psi-blast which maxmixe Sw measure on training set
./espritz.pl example_fastas/ pD 1
#execute the following for long disprot type disorder predictions using psi-blast at a 5% FPR
./espritz.pl example_fastas/ pD 0
#execute the following for NMR disorder predictions using psi-blast which maxmixe Sw measure on training set
./espritz.pl example_fastas/ pN 1
#execute the following for NMR disorder predictions using psi-blastusing psi-blast at a 5% FPR
./espritz.pl example_fastas/ pN 0

