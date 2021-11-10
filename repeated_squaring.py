from pyblq import Blq, ket

repeated_squaring = Blq("""[out:N @ x:2@@q <- x:2@@q ]

declare idx : q

repeat q:
    if x[idx]:
        declare tmp : N
        tmp += a
        repeat idx:
            tmp = tmp*tmp
        out += tmp
        discard tmp
    idx += 1

discard idx

""",a=7,N=15,q=8)



print(repeated_squaring * (ket(0).tensorpower(8))  )
