
def assert_equations(given, expected, diff=False):
    given_lines = sorted(given.__str__().replace(' ', '').strip().split('\n'))
    expected_lines = sorted(expected.__str__().replace(' ', '').strip().split('\n'))
    if diff:
        print("diff:")
        for g, e in zip(given_lines, expected_lines):
            print('given:    {}'.format(g))
            print('expected: {}'.format(e))
    assert given_lines == expected_lines
