#############################################################################
#############################################################################
#
# This file implements miscellaneous utility functions
#                                                                               
#############################################################################
#############################################################################


""" 
A generic error message for functionality that is currently unimplemented.

Useful for development and setting up type hierarchies.

This provides a fallback method if you forget to implement a function - e.g., if Julia determines
that it cannot find the relevant function for your Polynomial subtype, it will fallback to calling
a function that throws an error.

x:
    The input type
method: 
    The name of the unimplemented method
"""
function not_implemented_error(x, methodName::String)
    error("The method '$(methodName)' is not yet implemented for an object of type $(typeof(x))")
end