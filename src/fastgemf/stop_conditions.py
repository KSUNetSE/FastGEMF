def stop_cond(setup):
    """
    Determines if the simulation should stop based on the user-defined stopping condition.

    Parameters:
    -----------
    setup : SimulationSetu

    Returns:
    --------
    bool
        True if the stop condition is met, otherwise False.
    """

    if setup.stop_condition.get("events"):
        return (setup.iteration>setup.stop_condition.get("events") and len(setup.times) )
    
    elif setup.stop_condition.get("time"):
        return ( setup.Tf>setup.stop_condition.get("time") and len(setup.times) )
    
    else:
        return False
    