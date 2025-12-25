#!/usr/bin/ruby

N = 160
E = (N*0.8).to_i

puts "#pre,post_i,post_c,weight,tau,erev,delay,e/i"

N.times{|i|
  N.times{|j|
    weight = if i < E then 0.1e-3 else 0.4e-3 end
    tau = if i < E then 5 else 10 end
    erev = if i < E then 0.0 else -70.0 end
    ei = if i < E then 'e' else 'i' end
    puts "#{i},#{j},#{0},#{weight},#{tau},#{erev},1,#{ei}" if rand < 0.02
  }
}
