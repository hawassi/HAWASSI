function IQn = fun_IntegQn(kn,eta,D,z1,z2)
IQn=(sin(kn.*(D+z2))-sin(kn.*(D+z1)))./kn./cos(kn.*(D+eta));