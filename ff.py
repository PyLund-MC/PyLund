import Hepmcconverter
import Fourvector
import Particle

fv = Fourvector.vector4D(1,1,1,1)

ID = 2

part = Particle.particle4D(ID,fv)

print part.Getvector4D()