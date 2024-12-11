module InstantShell

struct Element 


end


struct Node 

end


struct Material{MaterialPropertyType<:Real}
    E_xx            ::MaterialPropertyType
    E_yy            ::MaterialPropertyType
    v_xx            ::MaterialPropertyType
    v_yy            ::MaterialPropertyType
    G_xy            ::MaterialPropertyType

    function Material(E_xx::Real, E_yy::Real, v_xx::Real, v_yy::Real, G_xy::Real)
        material_properties = promote(E_xx, E_yy, v_xx, v_yy, G_xy)

        return new{eltype(material_properties)}(material_properties...)
    end
end

#make a list, to identify all the elements at the beginning
#assign to the list





#assign material
#assign element type




end # module InstantShell
