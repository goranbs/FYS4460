#ifndef UNITCONVERTER_H
#define UNITCONVERTER_H

class UnitConverter{

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
    double to_aangstrom(double aangstrom);
    double from_aangstrom(double aangstrom);


private :
    static const double m_mass = 39.948; // amu
    static const double m_length = 3.405e-10; // m
    static const double m_time = 2.1569e-12; // s
    static const double m_energy = 1.0318*1.602e-20; // E
    static const double m_temperature = 119.74; // K
    static const double m_aangstrom = 3.405; // Å

};

#endif // UNITCONVERTER_H
