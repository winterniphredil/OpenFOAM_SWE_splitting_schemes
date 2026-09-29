class numericalParameters
{
// Public data
public:
    //const doubleScalar tolerance;
    const int maxIts;
    const Switch FB;
    
// Construct from dictionary
    numericalParameters(dictionary dict)
    :
    //tolerance(dict.lookup("tolerance")),
    maxIts(readLabel(dict.lookup("maxIts"))),
    FB(dict.lookup("forwardBackwardGravity"))
    {};
};

numericalParameters num(mesh.solution());
