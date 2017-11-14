//Nicolas Stoian
//CS780 Image Processing
//Project 9.1 - K-Curvature Corner Detector

//This program needs 5 command line arguments
//argv[1] "input1" for text file representing the boundary points of an object in an image
//argv[2] "input2" for K, the length of neighborhood in K-curvature computation
//argv[3] "output1" for The result of the K-curvature, as a text file representing the boundary points of the object with corner indicator.
//argv[4] "output2" for pretty print (displaying) the result of argv[3] as in an image.
//argv[5] "output3" for all debugging output

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace std;

class BoundryPt{ friend class KCurvature;
private:
    int x;
    int y;
    double curvature;
    int localMax;
    int corner;

public:
    BoundryPt();
    BoundryPt(int x, int y);
    ~BoundryPt();
    int getX();
    int getY();
    int getCorner();
};

class KCurvature{ friend class Image;;
private:
    int k;
    int numPts;
    BoundryPt** boundPtAry;
    int index;
    int q;
    int p;
    int r;

public:
    KCurvature();
    KCurvature(ifstream& inFile, int numPts, int k, ofstream& outFile);
    ~KCurvature();
    void loadData(int x, int y);
    double computeCurvature();
    void computeLocalMaxima();
    void isCorner();
    void printCurveInfo(ofstream& outFile);
    void printTextFile(ofstream& outFile, int numRows, int numCols, int minVal, int maxVal, int label);
};

class Image{
private:
    int** imgAry;
    int numRows;
    int numCols;
    int minVal;
    int maxVal;

public:
    Image();
    Image(KCurvature* imageCurvature, int numRows, int numCols, int minVal, int maxVal);
    ~Image();
    static string mapInt2Char(int theInt);
    void prettyPrint(ofstream& outFile);
};

int main(int argc, char* argv[]){
    ifstream inFile;
    inFile.open(argv[1]);
    int numRows, numCols, minVal, maxVal, label, numPts;
    inFile >> numRows;
    inFile >> numCols;
    inFile >> minVal;
    inFile >> maxVal;
    inFile >> label;
    inFile >> numPts;
    stringstream kValInput(argv[2]);
    int kValue;
    kValInput >> kValue;
    ofstream outFile3;
    outFile3.open(argv[5]);
    outFile3 << setw(2) << "Q" << "  " << setw(2) << "P" << "  " << setw(2) << "R" << "  " << setw(2)
             << "X" << "  " << setw(2) << "Y" << "  " << setw(9) << "Curvature" << endl << endl;
    KCurvature* imageCurvature = new KCurvature(inFile, numPts, kValue, outFile3);
    outFile3 << endl << endl << endl;
    imageCurvature->computeLocalMaxima();
    imageCurvature->isCorner();
    imageCurvature->printCurveInfo(outFile3);
    inFile.close();
    outFile3.close();
    ofstream outFile1;
    outFile1.open(argv[3]);
    imageCurvature->printTextFile(outFile1, numRows, numCols, minVal, maxVal, label);
    outFile1.close();
    Image* image = new Image(imageCurvature, numRows, numCols, minVal, maxVal);
    ofstream outFile2;
    outFile2.open(argv[4]);
    image->prettyPrint(outFile2);
    outFile2.close();
}

BoundryPt::BoundryPt(): x(0), y(0), curvature(0.0), localMax(0), corner(0){
}

BoundryPt::BoundryPt(int x, int y): curvature(0.0), localMax(0), corner(0){
    this->x = x;
    this->y = y;
}

BoundryPt::~BoundryPt(){
}

int BoundryPt::getX(){
    return x;
}

int BoundryPt::getY(){
    return y;
}

int BoundryPt::getCorner(){
    return corner;
}

KCurvature::KCurvature(): boundPtAry(NULL), k(0), numPts(0), index(0), q(0), p(0), r(0){
}

KCurvature::KCurvature(ifstream& inFile, int numPts, int k, ofstream& outFile){
    this->numPts = numPts;
    this->k = k;
    index = 0;
    boundPtAry = new BoundryPt*[numPts];
    int x, y;
    while(inFile >> x){
        inFile >> y;
        loadData(x, y);
        index++;
    }
    q = 0;
    p = k;
    r = 2 * k;
    do{
        index = p;
        double curvature = computeCurvature();
        boundPtAry[index]->curvature = curvature;
        outFile << setw(2) << q << "  " << setw(2) << p << "  " << setw(2) << r << "  " << setw(2) << boundPtAry[index]->x
                << "  " << setw(2) << boundPtAry[index]->y << "  " << setw(9) << boundPtAry[index]->curvature << endl;
        q = (q + 1) % numPts;
        p = (p + 1) % numPts;
        r = (r + 1) % numPts;
    }
    while(p != k);
}

KCurvature::~KCurvature(){
    if(boundPtAry != NULL){
        for(int i = 0; i < numPts; i++){
            delete boundPtAry[i];
        }
    }
    delete [] boundPtAry;
}

void KCurvature::loadData(int x, int y){
    boundPtAry[index] = new BoundryPt(x, y);
}

double KCurvature::computeCurvature(){
    double denom1 = ((double)boundPtAry[q]->x - (double)boundPtAry[p]->x);
    if(denom1 == 0.0){
        denom1 = denom1 + .00001;
    }
    double denom2 = ((double)boundPtAry[p]->x - (double)boundPtAry[r]->x);
    if(denom2 == 0.0){
        denom2 = denom2 + .00001;
    }
    return ((double)boundPtAry[q]->y - (double)boundPtAry[p]->y) / denom1 -
           ((double)boundPtAry[p]->y - (double)boundPtAry[r]->y) / denom2;
}

void KCurvature::computeLocalMaxima(){
    for(int i = 0; i < numPts; i++){
        if(fabs(boundPtAry[i]->curvature) >= fabs(boundPtAry[((i - 2) % numPts + numPts) % numPts]->curvature) &&
           fabs(boundPtAry[i]->curvature) >= fabs(boundPtAry[((i - 1) % numPts + numPts) % numPts]->curvature) &&
           fabs(boundPtAry[i]->curvature) >= fabs(boundPtAry[(i + 1) % numPts]->curvature) &&
           fabs(boundPtAry[i]->curvature) >= fabs(boundPtAry[(i + 2) % numPts]->curvature) &&
           boundPtAry[i]->curvature != 0){
            boundPtAry[i]->localMax = 1;
        }
    }
}

void KCurvature::isCorner(){
    for(int i = 0; i < numPts; i++){
        if(boundPtAry[i]->localMax){
            if(boundPtAry[((i - 1) % numPts + numPts) % numPts]->localMax &&
               boundPtAry[(i + 1) % numPts]->localMax){
                boundPtAry[i]->corner = 1;
            }
            else{
                boundPtAry[i]->corner = 8;
            }
        }
        else{
            boundPtAry[i]->corner = 1;
        }
    }
}

void KCurvature::printCurveInfo(ofstream& outFile){
    outFile << setw(5) << "Index" << "  " << setw(2) << "X" << "  " << setw(2) << "Y" << "  " << setw(9)
            << "Curvature" << "  " << setw(8) << "localMax" << "  " << setw(6) << "Corner" << endl << endl;
    for(int i = 0; i < numPts; i++){
        outFile << setw(5) << i << "  " << setw(2) << boundPtAry[i]->x << "  " << setw(2) << boundPtAry[i]->y << "  " << setw(9)
                << boundPtAry[i]->curvature << "  " << setw(8) << boundPtAry[i]->localMax << "  " << setw(6) << boundPtAry[i]->corner << endl;
    }
}

void KCurvature::printTextFile(ofstream& outFile, int numRows, int numCols, int minVal, int maxVal, int label){
    outFile << numRows << " " << numCols << " " << minVal << " " << maxVal << endl;
    outFile << label << endl;
    outFile << numPts << endl;
    for(int i = 0; i < numPts; i++){
        outFile << boundPtAry[i]->x << " " << boundPtAry[i]->y << " " << boundPtAry[i]->corner << endl;
    }
}

Image::Image(): imgAry(NULL), numRows(0), numCols(0), minVal(0), maxVal(0){
}

Image::Image(KCurvature* imageCurvature, int numRows, int numCols, int minVal, int maxVal){
    this->numRows = numCols; // These are reversed in the input file for a true (x,y) coordinate image representation, x = cols, y = rows.
    this->numCols = numRows; // I have reversed them here to get proper output but would recommend fixing the input file.
    this->minVal = minVal;
    this->maxVal = maxVal;
    imgAry = new int* [this->numRows];
    for(int i = 0; i < this->numRows; i++){
        imgAry[i] = new int [this->numCols];
    }
    for(int row = 0; row < this->numRows; row++){
        for(int col = 0; col < this->numCols; col++){
            imgAry[row][col] = 0;
        }
    }
    for(int i = 0; i < imageCurvature->numPts; i++){
        imgAry[imageCurvature->boundPtAry[i]->getY()][imageCurvature->boundPtAry[i]->getX()] = imageCurvature->boundPtAry[i]->getCorner();
    }
}

Image::~Image(){
    if(imgAry != NULL){
        for(int i = 0; i < numRows; i++){
            delete [] imgAry[i];
        }
    }
    delete [] imgAry;
}

string Image::mapInt2Char(int theInt){
    char toReturn [33];
    sprintf(toReturn, "%d", theInt);
    //itoa(theInt, toReturn, 10);
    return toReturn;
}

void Image::prettyPrint(ofstream& outFile){
    //for(int row = numRows - 1; row >= 0; row--){
    for(int row = 0; row < numRows; row++){
        for(int col = 0; col < numCols; col++){
            if(imgAry[row][col] <= 0){
                outFile << " " << " ";
            }
            else{
                outFile << mapInt2Char(imgAry[row][col]) << " ";
            }
        }
        outFile << endl;
    }
}
