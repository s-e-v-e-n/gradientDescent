#include <math.h>

const int inputs = 1;
const int tPoints = 3;
const int MIN = 0;
const int MID = 1;
const int MAX = 2;
const double step = 0.0007;
const double errorThreshold = 0.01;
const double iterations = 10;
const bool debug = false;
const bool debugcost = false;
const bool debughypothesis = false;
const bool debugfeaturescale = false;
const bool debugconvergence = false;

double (*valueMap)[tPoints] = new double[inputs][tPoints];   
double (*degreeMap)[tPoints] = new double[inputs][tPoints];   

void setup(){
  Serial.begin(9600);
}
void iout(bool debugflag, int value){
    if(debug==true||debugflag==true){
//        cout << value;
        Serial.print(value);
    }
}
void fout(bool debugflag, double value){
    if(debug==true||debugflag==true){
//        cout << value;
        Serial.print(value);
    }
}
void sout(bool debugflag, String value){
    if(debug==true||debugflag==true){
        if(value=="linebreak"){
        //          cout << "\n";
            Serial.println();
        }else{
        //          cout << value;
            Serial.print(value);
        }

    }
}

void loop(){
    Serial.println("Serialtest");
    sout(true,"string");
    sout(true,"linebreak");
    iout(true,1);
    sout(true,"linebreak");
    fout(true,1.5);
    sout(true,"linebreak");
    mainotbenamed();
    delay(1000);
}

void setMap(int input, int point, double value, double degree){
    valueMap[input][point] = value;
    degreeMap[input][point] = degree;
}

double featurescale(int input, double value, int power){
    double avg=0.0;
    double max=0.0;
    double min=10000000.0;
    double feat=0.0;
    // sum over all trainingsets
    for(int point=0; point<tPoints; point++){
        feat=pow(valueMap[input][point],power);
        if(feat > max){
            max=feat;
        }
        if(feat < min){
            min=feat;
        }
    }
    for(int point=0; point<tPoints; point++){
        feat=pow(valueMap[input][point],power);
        avg=avg+feat;
    }

    avg = avg/tPoints;
    fout(debugfeaturescale, value-avg);
    sout(debugfeaturescale, "\t ");
    fout(debugfeaturescale, value);
    sout(debugfeaturescale, "\t ");
    fout(debugfeaturescale, max-min);
    sout(debugfeaturescale, "linebreak");
    value = (value-avg)/(max-min);
    return value;
}

double hypothesis(int input, double value, double theta[4]){
    double xGuess = theta[0];
    double feat = 0.0;
    for(int t=1; t<4; t++){
        feat = pow(value,t);
        feat = featurescale(input,feat,t);
        xGuess = xGuess + (theta[t] * feat);
        sout(debughypothesis, " t: ");
        iout(debughypothesis, t);
        sout(debughypothesis, " feat: ");
        fout(debughypothesis, feat);
        sout(debughypothesis, "linebreak");
        sout(debughypothesis, " xGuess: ");
        fout(debughypothesis, xGuess);
        sout(debughypothesis, "linebreak");
    }
    return xGuess;
}

double costFunction(int input, double theta[4]){
    double sqrDiff = 0.0;
    double sum = 0.0;
    double xKnown = 0.0;                              // known value
    double xGuess = 0.0;                              // known value
    double yKnown = 0.0;
    double cost = 0.0;

    for(int point=0; point<tPoints; point++){ // known values
        xKnown = valueMap[input][point];
        yKnown = degreeMap[input][point];
        xGuess = hypothesis(input,xKnown,theta);
        sqrDiff = pow(xGuess - yKnown,2);
        sum = sum+sqrDiff;
    }
    cost = sum/2*tPoints;
    sout(debugcost,"cost:");
    fout(debugcost,cost);
    sout(debugcost,"\n");
    return cost;
}


double *gradientDescent(int input, double theta[4]){
    double thetaTemp[4];
    double cost = 1000000.0;
    double partDeriv;
    double oldcost = 0.0;
    double diff = 0.0;
    double sum = 0.0;
    double xKnown = 0.0;
    double xGuess = 0.0;
    double yKnown = 0.0;
    double xTra;                           

    for(int iteration=0; iteration<iterations; iteration++){
        // repeat until convergence
        oldcost = cost;
        cost=costFunction(input, theta);
        if(oldcost<cost){ 
            sout(debugconvergence,"-@iteration:" );
            iout(debugconvergence,iteration);
            sout(debugconvergence,"\n");
            break;
        }
        if(cost<errorThreshold){ 
            sout(debugconvergence,"t@iteration:" );
            iout(debugconvergence,iteration);
            sout(debugconvergence,"\n");
            break;
        }
        for(int t=1; t<4; t++){
            // sum over all trainingsets
            for(int point=0; point<tPoints; point++){
                xKnown = valueMap[input][point];
                yKnown = degreeMap[input][point];
                if(t > 0){
                    xTra = xKnown;
                }else{
                    xTra = 1;
                }
                xGuess = hypothesis(input,xKnown,theta);
                diff = (xGuess - yKnown) * xTra;
                sum = sum+diff;
            }
            partDeriv = sum/tPoints;
            thetaTemp[t] = theta[t] - step * partDeriv;
        }
        for(int t=0; t<4; t++){
            theta[t] = thetaTemp[t];
        }
    }
    return theta;
}

double getDegree(int input, double value){
    double degree =0.000;
    double theta[4];
    double *newtheta;
    theta[0]=0;
    theta[1]=0;
    theta[2]=0;
    theta[3]=0;
    newtheta = gradientDescent(input, theta);
    degree = hypothesis(input, value, newtheta);
    return degree;
}

int mainotbenamed (){
    double normalizedDegree;
    double normalizedValue;
    double offset;
    setMap(0,MIN,0.001,0.001);
    setMap(0,MID,0.900,0.090);
    setMap(0,MAX,1.024,0.360);
    for(int value = 1.0; value <= 1024.0; value += 1.0){
        normalizedValue=value/1000.0;
        normalizedDegree=getDegree(0,normalizedValue)*1000.0;
        if(value==1 && normalizedDegree<0){
            offset=0-normalizedDegree;
        } 
        fout(true, value);
        sout(true, ";"); 
        fout(true, normalizedDegree+offset);
        sout(true, "linebreak");
    }
    return 0;
}
