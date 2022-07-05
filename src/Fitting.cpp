#include "Arduino.h"
#include "Fitting.h"
#include <math.h>

void Interpolation::sizeUpdate(int _array_size)
{
    _array_size_ = _array_size;
}

Interpolation::Interpolation(Data _data[], int _array_size)
{
    _array_size_ = _array_size;
    for (int i = 0; i < _array_size_; i++)
    {
        _data_[i].x = _data[i].x;
        _data_[i].y = _data[i].y;
    }
    linear_availability = false;
}

void Interpolation::interpolationUpdate(Data _data[])
{
    for (int i = 0; i < _array_size_; i++)
    {
        _data_[i].x = _data[i].x;
        _data_[i].y = _data[i].y;
    }
    linear_availability = false;
}

// Lagrange's formula
double Interpolation::lagrange(double _target)
{
    double result = 0;
    for (int i = 0; i < _array_size_; i++)
    {
        double term = _data_[i].y;
        for (int j = 0; j < _array_size_; j++)
        {
            if (j != i)
                term = term * (_target - _data_[j].x) / (_data_[i].x - _data_[j].x);
        }
        result += term;
    }
    return result;
}

double Interpolation::linear(double _target)
{
    double maxArrayVal, minArrayVal;
    double sorting_queue_x[_array_size_];
    double a, b, c, d;
    for (int l = 0; l < _array_size_; l++)
    {
        sorting_queue_x[l] = _data_[l].x;
    }
    int key, j;
    for (int i = 1; i < _array_size_; i++)
    {
        key = sorting_queue_x[i];
        j = i - 1;
        while (j >= 0 && sorting_queue_x[j] > key)
        {
            sorting_queue_x[j + 1] = sorting_queue_x[j];
            j = j - 1;
        }
        sorting_queue_x[j + 1] = key;
    }
    maxArrayVal = sorting_queue_x[_array_size_ - 1];
    minArrayVal = sorting_queue_x[0];
    if ((_target > minArrayVal) && (_target < maxArrayVal))
    {
        for (int k = 0; k < _array_size_ - 1; k++)
        {
            if ((_target > sorting_queue_x[k]) && (_target < sorting_queue_x[k + 1]))
            {
                Index.x1 = k;
                Index.x2 = k + 1;
            }
        }
        for (int m = 0; m < _array_size_; m++)
        {
            if (_data_[m].x == sorting_queue_x[Index.x1])
            {
                Index.y1 = m;
            }
            if (_data_[m].x == sorting_queue_x[Index.x2])
            {
                Index.y2 = m;
            }
        }

        return ((_data_[Index.y2].y - _data_[Index.y1].y) / (sorting_queue_x[Index.x2] - sorting_queue_x[Index.x1])) * (_target - sorting_queue_x[Index.x2]) + _data_[Index.y2].y;
        linear_availability = DEFINED;
    }
    else
    {
        linear_availability = UNDEFINED;
    }
}

double Interpolation::getLinearDerivative(double _point)
{
    double maxArrayVal, minArrayVal;
    double sorting_queue_x[_array_size_];
    for (int l = 0; l < _array_size_; l++)
    {
        sorting_queue_x[l] = _data_[l].x;
    }
    int key, j;
    for (int i = 1; i < _array_size_; i++)
    {
        key = sorting_queue_x[i];
        j = i - 1;
        while (j >= 0 && sorting_queue_x[j] > key)
        {
            sorting_queue_x[j + 1] = sorting_queue_x[j];
            j = j - 1;
        }
        sorting_queue_x[j + 1] = key;
    }
    maxArrayVal = sorting_queue_x[_array_size_ - 1];
    minArrayVal = sorting_queue_x[0];
    if ((_point > minArrayVal) && (_point < maxArrayVal))
    {
        for (int k = 0; k < _array_size_ - 1; k++)
        {
            if ((_point > sorting_queue_x[k]) && (_point < sorting_queue_x[k + 1]))
            {
                Index.x1 = k;
                Index.x2 = k + 1;
            }
        }
        for (int m = 0; m < _array_size_; m++)
        {
            if (_data_[m].x == sorting_queue_x[Index.x1])
            {
                Index.y1 = m;
            }
            if (_data_[m].x == sorting_queue_x[Index.x2])
            {
                Index.y2 = m;
            }
        }
        return (_data_[Index.y2].y - _data_[Index.y1].y) / (sorting_queue_x[Index.x2] - sorting_queue_x[Index.x1]);
        linear_derivative_availability = DEFINED;
    }
    else
    {
        linear_derivative_availability = UNDEFINED;
    }
}

bool Interpolation::linearAvailable()
{
    return linear_availability;
}

bool Interpolation::linearDerivativeAvailable()
{
    return linear_derivative_availability;
}

double Interpolation::linearRegression(double _point)
{
    LinearRegression.sum_x1 = 0;
    LinearRegression.sum_x2 = 0;
    LinearRegression.sum_x1 = 0;
    LinearRegression.sum_y1 = 0;
    LinearRegression.sum_xy = 0;
    for (int i = 0; i < _array_size_; i++)
    {
        LinearRegression.sum_x1 = LinearRegression.sum_x1 + _data_[i].x;
        LinearRegression.sum_x2 = LinearRegression.sum_x2 + _data_[i].x * _data_[i].x;
        LinearRegression.sum_y1 = LinearRegression.sum_y1 + _data_[i].y;
        LinearRegression.sum_xy = LinearRegression.sum_xy + _data_[i].x * _data_[i].y;
    }

    LinearRegression.coeff = (_array_size_ * LinearRegression.sum_xy - LinearRegression.sum_x1 * LinearRegression.sum_y1) / (_array_size_ * LinearRegression.sum_x2 - LinearRegression.sum_x1 * LinearRegression.sum_x1);
    LinearRegression.constant = (LinearRegression.sum_y1 - LinearRegression.coeff * LinearRegression.sum_x1) / _array_size_;

    return _point * LinearRegression.coeff + LinearRegression.constant;
}

double Interpolation::getLinearRegressionDerivative()
{
    LinearRegression.sum_x1 = 0;
    LinearRegression.sum_x2 = 0;
    LinearRegression.sum_y1 = 0;
    LinearRegression.sum_xy = 0;
    for (int i = 0; i < _array_size_; i++)
    {
        LinearRegression.sum_x1 = LinearRegression.sum_x1 + _data_[i].x;
        LinearRegression.sum_x2 = LinearRegression.sum_x2 + _data_[i].x * _data_[i].x;
        LinearRegression.sum_y1 = LinearRegression.sum_y1 + _data_[i].y;
        LinearRegression.sum_xy = LinearRegression.sum_xy + _data_[i].x * _data_[i].y;
    }

    return (_array_size_ * LinearRegression.sum_xy - LinearRegression.sum_x1 * LinearRegression.sum_y1) / (_array_size_ * LinearRegression.sum_x2 - LinearRegression.sum_x1 * LinearRegression.sum_x1);
}

double Interpolation::exponentialFitting(double _point)
{
    Exponential.sum_x1 = 0;
    Exponential.sum_x2 = 0;
    Exponential.sum_y1 = 0;
    Exponential.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Exponential.sum_x1 = Exponential.sum_x1 + _data_[i].x;
        Exponential.sum_x2 = Exponential.sum_x2 + _data_[i].x * _data_[i].x;
        Exponential.sum_y1 = Exponential.sum_y1 + log(_data_[i].y);
        Exponential.sum_xy = Exponential.sum_xy + _data_[i].x * log(_data_[i].y);
    }

    Exponential.B = (_array_size_ * Exponential.sum_xy - Exponential.sum_x1 * Exponential.sum_y1) / (_array_size_ * Exponential.sum_x2 - Exponential.sum_x1 * Exponential.sum_x1);
    Exponential.A = (Exponential.sum_y1 - Exponential.B * Exponential.sum_x1) / _array_size_;

    Exponential.coeff = exp(Exponential.A);
    Exponential.base = exp(Exponential.B);

    return exp(Exponential.A) * pow(exp(Exponential.B), _point);
}

double Interpolation::getExponentialBase()
{
    Exponential.sum_x1 = 0;
    Exponential.sum_x2 = 0;
    Exponential.sum_y1 = 0;
    Exponential.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Exponential.sum_x1 = Exponential.sum_x1 + _data_[i].x;
        Exponential.sum_x2 = Exponential.sum_x2 + _data_[i].x * _data_[i].x;
        Exponential.sum_y1 = Exponential.sum_y1 + log(_data_[i].y);
        Exponential.sum_xy = Exponential.sum_xy + _data_[i].x * log(_data_[i].y);
    }

    Exponential.B = (_array_size_ * Exponential.sum_xy - Exponential.sum_x1 * Exponential.sum_y1) / (_array_size_ * Exponential.sum_x2 - Exponential.sum_x1 * Exponential.sum_x1);

    return exp(Exponential.B);
}

double Interpolation::getExponentialCoeff()
{
    Exponential.sum_x1 = 0;
    Exponential.sum_x2 = 0;
    Exponential.sum_y1 = 0;
    Exponential.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Exponential.sum_x1 = Exponential.sum_x1 + _data_[i].x;
        Exponential.sum_x2 = Exponential.sum_x2 + _data_[i].x * _data_[i].x;
        Exponential.sum_y1 = Exponential.sum_y1 + log(_data_[i].y);
        Exponential.sum_xy = Exponential.sum_xy + _data_[i].x * log(_data_[i].y);
    }

    Exponential.B = (_array_size_ * Exponential.sum_xy - Exponential.sum_x1 * Exponential.sum_y1) / (_array_size_ * Exponential.sum_x2 - Exponential.sum_x1 * Exponential.sum_x1);
    Exponential.A = (Exponential.sum_y1 - Exponential.B * Exponential.sum_x1) / _array_size_;

    return exp(Exponential.A);
}

double Interpolation::powerFitting(double _point)
{
    Power.sum_x1 = 0;
    Power.sum_x2 = 0;
    Power.sum_y1 = 0;
    Power.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Power.sum_x1 = Power.sum_x1 + log(_data_[i].x);
        Power.sum_x2 = Power.sum_x2 + log(_data_[i].x) * log(_data_[i].x);
        Power.sum_y1 = Power.sum_y1 + log(_data_[i].y);
        Power.sum_xy = Power.sum_xy + log(_data_[i].x) * log(_data_[i].y);
    }

    Power.power = (_array_size_ * Power.sum_xy - Power.sum_x1 * Power.sum_y1) / (_array_size_ * Power.sum_x2 - Power.sum_x1 * Power.sum_x1);
    Power.A = (Power.sum_y1 - Power.power * Power.sum_x1) / _array_size_;

    Power.coeff = exp(Power.A);

    return Power.coeff * pow(_point, Power.power);
}

double Interpolation::getPowerPower()
{
    Power.sum_x1 = 0;
    Power.sum_x2 = 0;
    Power.sum_y1 = 0;
    Power.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Power.sum_x1 = Power.sum_x1 + log(_data_[i].x);
        Power.sum_x2 = Power.sum_x2 + log(_data_[i].x) * log(_data_[i].x);
        Power.sum_y1 = Power.sum_y1 + log(_data_[i].y);
        Power.sum_xy = Power.sum_xy + log(_data_[i].x) * log(_data_[i].y);
    }

    return (_array_size_ * Power.sum_xy - Power.sum_x1 * Power.sum_y1) / (_array_size_ * Power.sum_x2 - Power.sum_x1 * Power.sum_x1);
}

double Interpolation::getPowerCoeff()
{
    Power.sum_x1 = 0;
    Power.sum_x2 = 0;
    Power.sum_y1 = 0;
    Power.sum_xy = 0;

    for (int i = 0; i < _array_size_; i++)
    {
        Power.sum_x1 = Power.sum_x1 + log(_data_[i].x);
        Power.sum_x2 = Power.sum_x2 + log(_data_[i].x) * log(_data_[i].x);
        Power.sum_y1 = Power.sum_y1 + log(_data_[i].y);
        Power.sum_xy = Power.sum_xy + log(_data_[i].x) * log(_data_[i].y);
    }

    Power.power = (_array_size_ * Power.sum_xy - Power.sum_x1 * Power.sum_y1) / (_array_size_ * Power.sum_x2 - Power.sum_x1 * Power.sum_x1);

    return exp((Power.sum_y1 - Power.power * Power.sum_x1) / _array_size_);
}