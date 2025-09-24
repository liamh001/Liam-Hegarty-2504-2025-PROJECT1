function test_z_mod_p()
    println("=== TESTING ZModP OPERATIONS ===\n")
    Random.seed!(0)
    
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    test_count = 0
    
    for i in 1:1000
        p = primes[rand(1:10)]
        a_val = rand(-10^4:10^4)
        b_val = rand(-10^4:10^4)
        
        a = ZModP(a_val, Val(p))
        b = ZModP(b_val, Val(p))
        
        @assert Int(a + b) == mod(a_val + b_val, p)
        @assert Int(a - b) == mod(a_val - b_val, p)
        @assert Int(a * b) == mod(a_val * b_val, p)
        @assert Int(-a) == mod(-a_val, p)
        
        if b != 0
            c = a ÷ b
            @assert Int(c * b) == Int(a)
        end
        
        n = rand(0:5)
        @assert Int(a^n) == mod(powermod(a_val, n, p), p)
        
        test_count += 1
        if i % 100 == 0
            print(".")
        end
    end
    
    println("\n✓ All $test_count ZModP tests passed")
    println("✓ Tested with primes: $primes")
    println("\n=== ALL ZModP TESTS PASSED ===")
end