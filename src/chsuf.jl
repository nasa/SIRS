"""
    chsuf(str, old_suffix, new_suffix)

Change a string's suffix

    Parameters: str::String
                  The string to operate on
                old_suffix::String
                  The existing suffix
                new_suffix::String
                  The desired new suffix
"""
function chsuf(str::String, old_suffix::String, new_suffix::String)
   return(str[1:findlast(old_suffix, str)[1]-1] * new_suffix)
end