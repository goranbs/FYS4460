#include "unitconverter.h"

UnitConverter::UnitConverter(){

}

double UnitConverter::to_amu(double mass){
    return m_mass*mass;
}
double UnitConverter::from_amu(double mass){
    return mass/m_mass;
}
double UnitConverter::to_meter(double length){
    return length*m_length;
}
double UnitConverter::from_meter(double length){
    return length/m_length;
}
double UnitConverter::to_time(double time){
    return time*m_time;
}
double UnitConverter::from_time(double time){
    return time/m_time;
}
double UnitConverter::to_energy(double energy){
    return energy*m_energy;
}
double UnitConverter::from_energy(double energy){
    return energy/m_energy;
}
