# for running on heroku redis

import os

import redis
from rq import Worker, Queue, Connection

# for redis server

listen = ['high', 'default', 'low']

redis_url = os.getenv('REDIS_URL', 'redis://:p2bc366791c653b891068ce502577528e16e4d711bec7ff6a2ff1fede43d2abd5@ec2-3-231-65-84.compute-1.amazonaws.com:20529')

conn = redis.from_url(redis_url)

if __name__ == '__main__':
    with Connection(conn):
        worker = Worker(map(Queue, listen))
        worker.work()


# for running on local machine

# import os

# import redis
# from rq import Worker, Queue, Connection

# listen = ['default']

# redis_url = os.getenv('REDIS_URL', 'redis://localhost:6379')

# conn = redis.from_url(redis_url)

# if __name__ == '__main__':
#     with Connection(conn):
#         worker = Worker(list(map(Queue, listen)))
#         worker.work()