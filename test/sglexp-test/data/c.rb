#!/usr/bin/ruby

N = 200
E = (N*0.8).to_i

puts "# pre,post_i,post_c,weight,tau,erev,delay,e/i"

N.times{|i|
  N.times{|j|
    weight = if i < E then 1.0 else 4.0 end
    tau = if i < E then 5.0 else 10.0 end
    erev = if i < E then 0.0 else -70.0 end
    ei = if i < E then 'e' else 'i' end
    puts "#{i},#{j},#{0},#{weight},#{tau},#{erev},1,#{ei}" if rand < 0.05
  }
}
