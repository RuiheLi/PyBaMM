#!/usr/bin/python3
import time, random
from multiprocessing import Pool   

# 这个类就是个筐，各种参数都可以往里装
class TaskConfig:
    def __init__(self, id, seconds):
        self.id = id
        # 运行参数
        self.seconds = seconds

        
# 运行结果（这里只考虑正常结果，也可以加上异常结果，用一个bool变量来区分成功和失败）
class TaskResult:
    def __init__(self, id, result):
        self.id = id
        self.result = result

        
# 模拟耗时很长的任务
def run_sim(config):
    time.sleep(config.seconds)
    result = 'Task %d finished after %d seconds' % (config.id, config.seconds) 
    return TaskResult(config.id, result)


# 生成配置的函数，每个任务随机取参数1~10   
def gen_config(id):
    return TaskConfig(id, random.randint(1,10))

    
def run_a_round(configs):
    worker_pool = Pool(processes=len(configs))  
    async_results = []
    for config in configs:
        async_result = worker_pool.apply_async(run_sim, (config,))
        async_results.append(async_result)
    # 不再允许提交任务
    worker_pool.close()
    # 设置5秒超时，主进程负责每1秒检查一次任务是否完成，注意检查间隔应该远小于超时时间，
    # 哪怕jupyter的sleep不准确也不会导致严重超时那种
    TIMEOUT = 5
    CHECK_INTERVAL = 1
    begin_time = time.time()
    while True:
        time.sleep(CHECK_INTERVAL)
        # 都完成了，就结束
        all_ready = all([async_result.ready() for async_result in async_results])
        if all_ready:
            break
        # 超时了，也结束
        current_time = time.time()
        if current_time - begin_time > TIMEOUT:
            break
            
    # 要么全部完成，要么超时
    task_results = []
    for i, async_result in enumerate(async_results):
        if async_result.ready():
            task_result = async_result.get()
            print(task_result.result)
            task_results.append(task_result)
        else:
            print('Task %d timeout' % i)
    
    # 关闭进程池释放资源
    worker_pool.terminate()
    
    return task_results

    
if __name__ == '__main__': 
    # 初始化10个任务
    task_configs = [gen_config(i) for i in range(10)]
    
    # 跑到没有任务可跑为止，最多跑三轮
    round_cnt = 1
    while len(task_configs) > 0 and round_cnt <= 3:
        print('--- round %d ---' % round_cnt)
        task_results = run_a_round(task_configs)
        # 根据结果更新下一轮参数
        task_configs = [gen_config(task_result.id) for task_result in task_results]
        
        round_cnt += 1
    
    
    