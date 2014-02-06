#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H

class UnitConverter
{
public:
    UnitConverter();
    double to_amu(double mass);
    double from_amu(double mass);
    double to_meter(double length);
    double from_meter(double length);
    double to_time(double time);
    double from_time(double time);
    double to_energy(double energy);
    double from_energy(double energy);
    }

private :
    const double m_mass = 39.948; // amu
    const double m_length = 3.405e-10; // m
    const double m_time = 2.1569e-12; // s
    const double m_energy = 1.0318*1.602e-20; // E
    const double m_temperature = 119.74; // K
};

#endif // UNITCONVERTER_H
