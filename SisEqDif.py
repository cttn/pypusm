def PuSm_unidades(concentraciones, t, flujo, constantes):
    """ Sistema de ecuaciones para el transitorio de Plutonio-Samario. Las 
    concentraciones de nucleidos y el flujo neutrónico deben darse en 
    unidades 'absolutas'"""

    Csm, Cpm, Cnp, Cpu = concentraciones
    Ssm, Spu, Lpm, Lnp = constantes

    eqs =  [Lpm*Cpm - Ssm*flujo*Csm,
           -Lpm*Cpm + Lpm*flujo,
           -Lnp*Cnp + Lnp*flujo,
            Lnp*Cnp - Spu*flujo*Cpu]
    return eqs

def pusm_re(concentraciones, t, FPP):
    """ Sistema de ecuaciones para el transitorio de Plutonio-Samario. Las 
    concentraciones de nucleídos deben darse como fracción respecto de las
    concentraciones en el flujo de equilibrio, 1FPP. El flujo es una función 
    del tiempo, como fraccion de plena potencia"""

    #Esm, Epm, Enp, Epu = [1.0, 1.0, 1.0, 1.0] #concentraciones INICIALES de equilibrio
    Esm, Epm, Enp, Epu = concentraciones
    Ssm, Spu, Lpm, Lnp = [2.0*4.1E-20, 1.5*1.0113E-21, 0.313416, 0.294336] #[cm^2,cm^2,d^-1,d^-1]
    flujo_eq = 86400*1.0E+14 # 86400s/d*n/s·cm^2 = n/d·cm^2
    
    eqs = [ Ssm * flujo_eq * ( Epm - FPP(t) * Esm ),
            Lpm * ( FPP(t) - Epm ),
            Lnp * ( FPP(t) - Enp ),
            Spu * flujo_eq * ( Enp - FPP(t) * Epu )]
    return eqs



