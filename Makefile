DATES = S19980012010031 S19980322010059 S19980602010090 S19980912010120 \
        S19981212010151 S19981522010181 S19981822010212 S19982132010243 \
        S19972442010273 S19972742010304 S19973052010334 S19973352010365
FILES = $(foreach d,$(DATES),$(d).L3m_MC_CHL_chlor_a_9km.nc)
TOPO = GEBCO_08_v1.nc
TARG = seawifs-clim-1997-2010.nc
HASH = hash.md5

# Main target checks the checksums
check: $(HASH)
	md5sum -c $(HASH)
Check: check
	md5sum -c raw.md5
	md5sum -c gebco.md5

# Fetch data from web
S%.nc:
	wget https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/$@

# Fetch topography file from archive (not easily available on web)
GEBCO_08_v1.nc:
	cp /archive/gold/datasets/topography/$@ .
	md5sum -c gebco.md5
seawifs_ocn_mask.nc:
	wget -nv -O $@.gz ftp://ftp.gfdl.noaa.gov/pub/aja/seawifs_ocn_mask.nc.gz
	gunzip $@.gz
	md5sum -c seawifs_ocn_mask.nc.md5

# Fill and join data
$(TARG): $(TOPO) $(FILES)
	./fill_and_join_chlor_a.py $@ $(TOPO) $(FILES)

stage.%.nc: %.L3m_MC_CHL_chlor_a_9km.nc
	./fill_and_join_chlor_a.py $@ $(TOPO) $^
stages: $(foreach d,$(DATES),stage.$(d).nc)
	ncrcat -h $^ test.nc
test: stage.S19980012010031.nc
	md5sum $^

# Record checksums
$(HASH): | $(TARG)
	md5sum $(TARG) > $@
raw.md5: $(FILES)
	md5sum $(FILES) > $@

# Clean up
clean:
	rm *.nc
