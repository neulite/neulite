#!/usr/bin/ruby

N = 160
E = (N*0.8).to_i

puts "#pre,post_i,post_c,weight,decay,rise,erev,delay,e/i"

N.times{|i|
  N.times{|j|
    weight = if i < E then 0.1e-3 else 0.4e-3 end
    tau_decay = if i < E then 5 else 10 end
    tau_rise = 1
    erev = if i < E then 0.0 else -70.0 end
    ei = if i < E then 'e' else 'i' end
    puts "#{i},#{j},#{0},#{weight},#{tau_decay},#{tau_rise},#{erev},1,#{ei}" if rand < 0.02
  }
}
