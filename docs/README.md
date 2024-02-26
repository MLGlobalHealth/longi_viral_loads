# Data Preprocessing

The flowchart below illustrates in detail how the data is consumed and processed by various R scripts at each stage. The coloured silos, rounded boxes, and grey boxes, represent data, R scripts, and directories respectively. The data is coloured according to the [medallion architecture](https://www.databricks.com/glossary/medallion-architecture), where bronze represents raw data, silver represents cleansed and processed data, and gold represents data which are consumption ready (i.e. directly used within statistical models).

```mermaid
flowchart LR 
	subgraph data
        d11[(path.flow.r1518)]
        d12[(path.flow.r19)]
        d23[(path.negatives.r1520)]
        d24[(path.nm.participants)]
        d25[(path.participation)]
        d26[(path.viral.loads1)]
        d27[(path.viral.loads2)]

	end 
	style d11 fill:#CD7F32
	style d12 fill:#CD7F32
	style d23 fill:#CD7F32
	style d24 fill:#CD7F32
	style d25 fill:#CD7F32
	style d26 fill:#CD7F32
	style d27 fill:#CD7F32
	
    s1(src/01_get_census_eligible_count.R)
    s2(src/02_new_process_data.R)
    s3(src/03_get_participation_rates.R) 
    s4(scripts/VL_run_cmdstan.R)

    subgraph Processed Data
        subgraph Confidential Data
            c11[(path.hivstatusvl.r1520)]
            style c11 fill:#F5F5F5
            style c11 color:#000000
        end
        subgraph Zenodo data
            o11[(path.comm.censsize)]
            o12[(path.census.eligible.aggregated)]
            o21[(path.participation.rates)]
            o31[(path.aggregated.nums.denoms.r1619)]

            style o11 fill:#F5F5F5
            style o12 fill:#F5F5F5
            style o21 fill:#F5F5F5
            style o31 fill:#F5F5F5
            style o11 color:#000000
            style o12 color:#000000
            style o21 color:#000000
            style o31 color:#000000
        end
    end
	

	
	d11 -- load --> s1
	d12 -- load --> s1
	s1 -- create --> o11
	s1 -- create --> o12
  
	d11 -- load --> s2
	d12 -- load --> s2
	d23 -- load --> s2
	d24 -- load --> s2
	d25 -- load --> s2
	d26 -- load --> s2
	d27 -- load --> s2
	s2 -- create --> c11
 
	c11 -- load --> s3
	o12 -- load --> s3
	s3 -- create --> o21
    
    c11 -- load --> s4
    s4 -- create --> o31
```


Successively, the outputs of the preprocessing are used in order to run our models.

```mermaid
flowchart LR 
	subgraph Zenodo data
		o11[(path.aggregated.nums.denoms.r1619)]
        o21[(path.participation.rates)]

    style o11 fill:#F5F5F5
    style o21 fill:#F5F5F5
    style o11 color:#000000
    style o21 color:#000000
	end 
 
	sc1(scripts/VL_run_stan.R)
    sc2(scripts/VL_postprocessing.R)
	
	subgraph results
		r11[(results/run-gp-*/cmd_*.rds)]
        r21[(results/tables/*)]

		style r11 fill:#ffd700
		style r21 fill:#ffd700
		style r11 color:#000000
		style r21 color:#000000
	end
	
	o11 -- load --> sc1
	o21 -- load --> sc1
	sc1 -- create --> r11

    r11 -- load --> sc2
    sc2 -- create --> r21
```
