using LinearAlgebra
using StaticArrays
using Cartesian
using Test

@testset "Cartesian.jl" begin
    @show @test Vector3()===SVector{3}([0.0, 0.0, 0.0])
    @show @test Vector3(1,2,3)===SVector{3}([1.0, 2.0, 3.0])
    @show @test Vector3([1,2,3])===SVector{3}([1.0, 2.0, 3.0])

    @show @test Matrix3()===SMatrix{3,3}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
    @show @test Matrix3(0*I)===SMatrix{3,3}([0.0 0.0 0.0; 0.0 0.0 0.0; 0.0 0.0 0.0])
    @show @test Matrix3(1*I)=== SMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    @show @test Matrix3(1,2,3,4,5,6,7,8,9) === SMatrix{3,3}([1.0 2.0 3.0; 4.0  5.0  6.0; 7.0  8.0  9.0])

    @show @test ô === Vector3()
    @show @test î === Vector3(1,0,0)
    @show @test ĵ === Vector3(0,1,0)
    @show @test k̂ === Vector3(0,0,1)

    @show @test Ô === Matrix3(0*I)
    @show @test Î === Matrix3(1*I)

    @show r = Vector3(1,2,3)
    @show @test r.X == 1 && r.Y == 2 && r.Z == 3
    @show @test r.SumSquares == 14
    @show @test r.Magnitude == sqrt(14)
    @show n = normalize(r)
    @show @test n == SVector{3}([1.0,2.0,3.0]/sqrt(14))
    @show @test n.Magnitude == 1.0

    @show g = Vector3(-2,3,2)

    @show @test dot(r,g)==10 && r⋅g == 10
    @show @test cross(r,g)===Vector3(-5,-8,7) && r×g === Vector3(-5,-8,7)

    @show @test cross(r)*g === cross(r,g)
    @show @test cross(r) === Matrix3(0,-3,2, 3,0,-1, -2,1,0) && ×(r) === cross(r)
    @show @test cross2(r) === cross(r)*cross(r) && ×(r)*×(r) === cross2(r)
    @show @test -×(r)*×(g) === Matrix3(12, 4, 6, -3, 4, -9, -2, -4, 4)

    @show A = Matrix3(5,-1,2,-1,5,3,-1,0,3)
    @show @test isapprox(85*(A\r) , Vector3(-18, -17, 79))
    @show @test isapprox( 85*inv(A), Matrix3(15,3,-13, 0, 17, -17, 5, 1, 24))
end
