#include <Fitting.h>
#include <Arduino.h>

Data func[] = {{0, 3}, {1, 6}, {2, 12}, {3, 24}, {4, 48}};
//Data k[];
const float example_data = 13.5;

Interpolation Test(func, sizeof(func) / sizeof(func[0]));

void setup()
{
    Serial.begin(9600);
    
    //Stream.sizeUpdate(array_size);
    //Stream.interpolationUpdate(k)
}

void loop()
{
    Serial.println(Test.exponentialFitting(example_data));
    Serial.println(Test.lagrange(example_data));
    Serial.println(Test.linearRegression(example_data));
    Serial.println(Test.powerFitting(example_data));
    Serial.println(Test.getExponentialBase());
}