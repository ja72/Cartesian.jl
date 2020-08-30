"""

Definitions for cartesian vectors and matrices used in 3D mechanics. 
The repository is located at https://github.com/ja72/Cartesian.jl

    $(EXPORTS)

    $(SIGNATURES)

    $(METHODLIST)

"""
module Cartesian
using LinearAlgebra
using StaticArrays
using Rotations
using DocStringExtensions

using Base: @propagate_inbounds

import Base.convert, Base.*
import LinearAlgebra.cross, LinearAlgebra.dot

export 
    Vector3,
    Matrix3,
    dot,
    cross,
    cross2,
    inv,
    solve

    """
    Cartesian 3×1 vector

        Vector3()
        Vector3(x,y,z)
        Vector3(vector)

    """
    const Vector3 = SVector{3, Float64}
    """
    Cartesian 3×3 matrix

        Matrix3()
        Matrix3(a11,a12,a13,a21,a23,a23,a31,a32,a33)
        Matrix3(matrix)

    """
    const Matrix3 = SMatrix{3, 3, Float64} 

    Vector3() = Vector3(0.0,0.0,0.0) 
    Matrix3() = Matrix3(0*I) 
    Matrix3(a11,a12,a13,a21,a22,a23,a31,a32,a33) = Matrix3([a11 a12 a13; a21 a22 a23; a31 a32 a33])

    function Base.getproperty(v::Vector3, name::Symbol)
        if name===:X  
            return v[1]
        end
        if name===:Y  
            return v[2]
        end
        if name===:Z
            return v[3]
        end
        if name===:SumSquares
            return dot(v,v)
        end
        if name===:Magnitude
            return sqrt(dot(v,v))
        end
        return getfield(v, name)
    end

    function Base.propertynames(v::Vector3)
        return (:X, :Y, :Z, :Magnitude, :SumSquares)
    end

    # function Base.getproperty(m::Matrix3, name::Symbol)
    #     if name===:A11; return m[1,1]; end
    #     if name===:A12; return m[1,2]; end
    #     if name===:A13; return m[1,3]; end
    #     if name===:A21; return m[2,1]; end
    #     if name===:A22; return m[2,2]; end
    #     if name===:A23; return m[2,3]; end
    #     if name===:A31; return m[3,1]; end
    #     if name===:A32; return m[3,2]; end
    #     if name===:A33; return m[3,3]; end
    # end

    # function Base.propertynames(m::Matrix3)
    #     return (:A11, :A12, :A13, :A21, :A22, :A23, :A31, :A32, :A33)
    # end

    export ô, î, ĵ, k̂, Ô, Î

    """
    Origin vector ``\\hat{o}=\\pmatrix{0\\\\0\\\\0}``

        const ô = Vector3(0,0,0)

    """
    const ô = Vector3(0,0,0)
    """
    X-axis basis vector ``\\hat{\\imath}=\\pmatrix{1\\\\0\\\\0}``

        const î = Vector3(1,0,0)

    """
    const î = Vector3(1,0,0)
    """
    Y-axis basis vector ``\\hat{\\jmath}=\\pmatrix{0\\\\1\\\\0}``

        const ĵ = Vector3(0,1,0)

    """
    const ĵ = Vector3(0,1,0)
    """
    Z-axis basis vector ``\\hat{k}=\\pmatrix{0\\\\0\\\\1}``

        const k̂ = Vector3(0,0,1)

    """
    const k̂ = Vector3(0,0,1)

    """
    Zero 3×3 matrix ``\\hat{O} = \\pmatrix{0 & 0 & 0 \\\\ 0 & 0 & 0 \\\\ 0 & 0 & 0}``

        const Ô = Matrix3(
            0,0,0,
            0,0,0,
            0,0,0)

    """
    const Ô = Matrix3()
    """
    Identity 3×3 matrix ``\\hat{I} = \\pmatrix{1 & 0 & 0 \\\\ 0 & 1 & 0 \\\\ 0 & 0 & 1}``

        const Î = Matrix3(
            1,0,0,
            0,1,0,
            0,0,1)

    """
    const Î = Matrix3(I)

    """
    Returns a normalized a vector whose magnitude is 1.0 or 0.0

        n = normalize(v)

    """
    function normalize(v::Vector3)::Vector3
        m2 = dot(v,v)
        if m2>0 && m2!=1.0
            return v/sqrt(m2)
        else
            return v
        end
    end

    """
    Dot product between vectors ``\\pmatrix{x \\\\ y \\\\ z} \\cdot \\pmatrix{u \\\\ v \\\\ w} = x*u + y*v + z*w``

        dot(a::Vector3, b::Vector3)
        a ⋅ b => dot(a,b)
        
    """
    @inline dot(a::Vector3, b::Vector3) = a[1]*b[1]+a[2]*b[2]+a[3]*b[3]

    """
    Cross product of vectors

        cross(a::Vector3, b::Vector3)::Vector3
        a × b => cross(a,b)

    """
    cross(a::Vector3, b::Vector3) = Vector3(
        a[2]*b[3]-a[3]*b[2],
        a[3]*b[1]-a[1]*b[3],
        a[1]*b[2]-a[2]*b[1])

    """
    Cross product matrix operator ``\\pmatrix{x \\\\ y \\\\ z}\\times =  \\pmatrix{0 & -z & y \\\\ z & 0 & -x \\\\ -y & x & 0}``

        cross(a::Vector3)
        ×(a::Vector3)

    """
    cross(a::Vector3) = Matrix3(0,-a[3],a[2],a[3],0,-a[1],-a[2],a[1],0)

    """
    Double cross product matrix operator ``\\pmatrix{x \\\\ y \\\\ z}\\times \\pmatrix{x \\\\ y \\\\ z}\\times =  \\pmatrix{-y^2-z^2 & x y & x z \\\\ x y & -x^2-z^2 & y z \\\\ x z & y z & -x^2-z^2}``

        cross2(a::Vector3)

    """
    cross2(a::Vector3) = Matrix3( 
        -a[2]^2-a[3]^2, a[1]*a[2], a[1]*a[3],  
        a[1]*a[2], -a[1]^2-a[3]^2, a[2]*a[3],
        a[1]*a[3], a[2]*a[3], -a[1]^2-a[2]^2)

    function Base.inv(m::Matrix3)::Matrix3
        t2 = m[1,1]*m[2,2]*m[3,3]
        t3 = m[1,2]*m[2,3]*m[3,1]
        t4 = m[1,3]*m[2,1]*m[3,2]
        t7 = m[1,1]*m[2,3]*m[3,2]
        t8 = m[1,2]*m[2,1]*m[3,3]
        t9 = m[1,3]*m[2,2]*m[3,1]
        t6 = 1.0/(t2+t3+t4-t7-t8-t9)
        return Matrix3(
             t6*(m[2,2]*m[3,3]-m[2,3]*m[3,2]),
            -t6*(m[1,2]*m[3,3]-m[1,3]*m[3,2]),
             t6*(m[1,2]*m[2,3]-m[1,3]*m[2,2]),
            -t6*(m[2,1]*m[3,3]-m[2,3]*m[3,1]),
             t6*(m[1,1]*m[3,3]-m[1,3]*m[3,1]),
            -t6*(m[1,1]*m[2,3]-m[1,3]*m[2,1]),
             t6*(m[2,1]*m[3,2]-m[2,2]*m[3,1]),
            -t6*(m[1,1]*m[3,2]-m[1,2]*m[3,1]),
             t6*(m[1,1]*m[2,2]-m[1,2]*m[2,1]))          
    end

    function solve(m::Matrix3, v::Vector3)::Vector3
        t = (-m[1,2]*m[3,3]+m[1,3]*m[3,2])*m[2,1]+(-m[2,2]*m[1,3]+m[2,3]*m[1,2])*m[3,1]+(m[2,2]*m[3,3]-m[2,3]*m[3,2])*m[1,1]
        return Vector3(
            (v[1]*( m[2,2]*m[3,3] - m[2,3]*m[3,2]) + v[2]*(-m[1,2]*m[3,3] + m[1,3]*m[3,2]) + v[3]*(-m[2,2]*m[1,3] + m[2,3]*m[1,2]))/t,
            (v[1]*(-m[2,1]*m[3,3] + m[2,3]*m[3,1]) + v[2]*( m[1,1]*m[3,3] - m[1,3]*m[3,1]) + v[3]*( m[2,1]*m[1,3] - m[2,3]*m[1,1]))/t,
            (v[1]*( m[2,1]*m[3,2] - m[2,2]*m[3,1]) + v[2]*(-m[1,1]*m[3,2] + m[1,2]*m[3,1]) + v[3]*(-m[2,1]*m[1,2] + m[2,2]*m[1,1]))/t) 
    end

    function Base.:\(A::Matrix3, b::Vector3)
        return solve(A,b)
    end    

end 
    
