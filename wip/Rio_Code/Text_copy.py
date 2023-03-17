                i_try = 0
                while i_try<3:
                    try:
                        getSth = Sol[key].entries[-1]
                    except:
                        i_try += 1
                        print(f"Fail to read Throughput capacity for the {i_try}th time")
                    else:
                        break