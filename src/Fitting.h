#ifndef Fitting_h
#define Fitting_h

#define UNDEFINED false
#define DEFINED true

struct Data
{
    double x, y;
};

struct 
{
    int x1, x2,y1,y2;
} Index;
 
struct 
{
    double sum_x1, sum_x2, sum_y1, sum_xy, constant, coeff;
} LinearRegression;

struct
{
    double sum_x1, sum_x2, sum_y1, sum_xy, coeff, base, A, B;
} Exponential;

struct
{
    double sum_x1, sum_x2, sum_y1, sum_xy, coeff, power, A;
} Power;

class Interpolation
{
    protected:
        //Availability of
        bool linear_availability = false;
        bool linear_derivative_availability = false;
        int _array_size_;

    public:
        Interpolation(Data _data[], int _array_size);
        void interpolationUpdate(Data _data[]);
        void sizeUpdate(int _array_size);

        double linear(double _target);
        double getLinearDerivative(double _point);
        bool linearDerivativeAvailable();
        bool linearAvailable();

        double lagrange(double _target);

        double linearRegression(double _point);
        double getLinearRegressionDerivative();

        double exponentialFitting(double _point);
        double getExponentialBase();
        double getExponentialCoeff();

        double powerFitting(double _point);
        double getPowerCoeff();
        double getPowerPower();

    private:
        Data _data_[];
};

#endif