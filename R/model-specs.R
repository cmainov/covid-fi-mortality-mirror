##  jhu models ##

# The Null Model
f.jhu.model.null <- deaths.jhu.adj ~ f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age ("Basic Model") ##
f.jhu.model.2 <- deaths.jhu.adj ~ I( fi.perc.20 / 4 ) + health.index + median.age + sir.jhu + f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# The "Basic Model" with some State-Level Variables, ("Basic + State Model")
f.jhu.model.3 <- deaths.jhu.adj ~ I( fi.perc.20 / 4 ) + health.index + median.age + sir.jhu +
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Add Other Demographic Variables ("The Demographic Model")
f.jhu.model.4 <- deaths.jhu.adj ~ I( fi.perc.20 / 4 ) + health.index + sir.jhu +
  median.age + perc.asian + perc.nh.white + perc.female + perc.native + 
  pop.density + disability + ed.1less.than.hspct + ed.5college.plus.pct +
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Add the rest of the County-Level Variables (The "Final Model")
f.jhu.model.5 <- deaths.jhu.adj ~ I(fi.perc.20/4) + pct.emp.trans + no.vehic + 
  disability + no.health.insur + perc.female + perc.nh.white +
  perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.jhu + urb.cat.code + 
  sir.jhu.state + health.index.state + elec.2020.margin.state +
  perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

f.jhu.model.5.pov <- deaths.jhu.adj ~ I(fi.perc.20/4) + pct.emp.trans + unemp.rate + no.vehic + 
  disability + no.health.insur + perc.black + perc.female + perc.hisp + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.jhu + urb.cat.code + 
  sir.jhu.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g, scale.model = T ) +
  f( state, model = "iid") + f( fips, model = "iid" )


##  cdc models ##

# The Null Model
f.cdc.model.null <- deaths.cdc ~ f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Models with Random Effects, Offset Terms, and Fixed Effects for Health Index and Median Age ("Basic Model") ##
f.cdc.model.2 <- deaths.cdc ~ I( fi.perc.20 / 4 ) + health.index + median.age + sir.cdc + f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# The "Basic Model" with some State-Level Variables, ("Basic + State Model")
f.cdc.model.3 <- deaths.cdc ~ I( fi.perc.20 / 4 ) + health.index + median.age + sir.cdc +
  sir.cdc.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Add County-level FI Fixed Effect and Other Demographic Variables ("The Demographic Model")
f.cdc.model.4 <- deaths.cdc ~ I( fi.perc.20 / 4 ) +
  median.age + perc.asian + perc.nh.white + perc.female + perc.native + 
  pop.density + disability + ed.1less.than.hspct + ed.5college.plus.pct +
  sir.cdc.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )

# Add the rest of the County-Level Variables (The "Final Model")
f.cdc.model.5 <- deaths.cdc ~ I(fi.perc.20/4) + pct.emp.trans + no.vehic + 
  disability + no.health.insur + perc.female + perc.nh.white +
  perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.cdc + urb.cat.code + 
  sir.cdc.state + health.index.state + elec.2020.margin.state +
  perc.vaccinated.state +
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )


f.cdc.model.5.pov <-deaths.cdc ~ I(fi.perc.20/4) + pct.emp.trans + unemp.rate + no.vehic + 
  disability + no.health.insur + perc.black + perc.female + perc.hisp + 
  perc.nh.white + perc.native + pop.density + ed.1less.than.hspct + 
  ed.5college.plus.pct + pct.emp.trade + median.age + perc.asian +
  poverty.rate + perc.vaccinated + gini.index + avg.hhsize + elec.2020.margin +
  ratio.pop.edp + health.index + sir.cdc + urb.cat.code + 
  sir.cdc.state + health.index.state + elec.2020.margin.state + perc.vaccinated.state +
  f( re.s, model = "besag", graph = g.cdc, scale.model = T ) +
  f( state, model = "iid" ) + f( fips, model = "iid" )