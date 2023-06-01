    bool isothermal = false;
    bool ENSE = true;

    mesh.solutionDict().subDict("NSE").readIfPresent
    (
        "isothermal",
        isothermal
    );

    mesh.solutionDict().subDict("NSE").readIfPresent
    (
        "ENSE",
        ENSE
    );
