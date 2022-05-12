import random

PopProfile_T = list[list[float]]

def VecDistance(a, b):
    return sum([abs(a[i] - b[i]) for i in range(0, len(a))])

def maxmin_full_info_sample(
    pop_performances: PopProfile_T,
    sample_size: int,
    init_test_id: int = -1
):
    assert len(pop_performances) > 0
    assert len(pop_performances[0]) > 0

    num_tests = len(pop_performances[0])
    assert init_test_id < num_tests

    test_profiles = [[pop_performances[org_id][test_id] for org_id in range(len(pop_performances))] for test_id in range(num_tests)]

    # If initial test not chosen, pick randomly.
    if init_test_id < 0:
        init_test_id = random.randrange(0, num_tests)

    available_tests = {test_id for test_id in range(num_tests)}
    included_tests = set()

    available_tests.remove(init_test_id)
    included_tests.add(init_test_id)

    # Sample until we've reached requested sample size.
    while len(included_tests) < sample_size:
        maxmin_dist = None
        maxmin_id = None

        # For each test we could sample, find its minimum distance to something already sampled.
        avail_order = list(available_tests)
        random.shuffle(avail_order)
        for avail_test_id in avail_order:
            cur_min_dist = None
            cur_test_profile = test_profiles[avail_test_id]

            for sampled_test_id in included_tests:
                sampled_test_profile = test_profiles[sampled_test_id]
                dist = VecDistance(cur_test_profile, sampled_test_profile)
                if (cur_min_dist==None):
                    cur_min_dist = dist
                elif (dist<cur_min_dist):
                    cur_min_dist = dist
                if cur_min_dist == 0: break

            if (maxmin_dist==None):
                maxmin_dist = cur_min_dist
                maxmin_id = avail_test_id
            elif (cur_min_dist > maxmin_dist):
                maxmin_dist = cur_min_dist
                maxmin_id = avail_test_id

        available_tests.remove(maxmin_id)
        included_tests.add(maxmin_id)

    return included_tests


pop_profile = [
    [1, 1, 0, 0, 0, 0],
    [1, 1, 1, 0, 0, 0],
    [1, 0, 0, 0, 1, 0],
    [1, 0, 0, 0, 0, 1],
]

result = maxmin_full_info_sample(
    pop_performances=pop_profile,
    sample_size=3
)

print(result)