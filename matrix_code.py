def det(a):
    if len(a) == 1:
        return a[0][0]
    if len(a) == 2:
        return a[0][0]*a[1][1] - a[0][1]*a[1][0]
    total = 0
    for i in range(len(a)):
        total += ((-1)**i)*a[0][i]*det([row[:i] + row[i+1:] for row in a[1:]])
    return total


def has_inverse(a):
    determinant = det(a)
    if abs(determinant) < 1e-10:
        return False
    return True


def inverse_with_det(a):
    determinant = det(a)
    if len(a) == 1:
        return [[1.0/a[0][0]]]
    if len(a) == 2:
        return [[a[1][1]/determinant, -a[0][1]/determinant],
                [-a[1][0]/determinant, a[0][0]/determinant]]
    cofactors = []
    for r in range(len(a)):
        cofRow = []
        for c in range(len(a)):
            minor  = [row[:c] + row[c+1:] for i, row in enumerate(a) if i != r]
            cofRow.append(((-1)**(r+c))*det(minor))
        cofactors.append(cofRow)
    adj = [[cofactors[c][r] for c in range(len(a))] for r in range(len(a))]
    
    return [[adj[r][c]/determinant for c in range(len(a))] for r in range(len(a))]


def gauss_jordan_inverse(a):
    p = [[0 for _ in range(len(a))]for _ in range(len(a))]
    for i in range(len(a)):
        p[i][i] = 1
    for i in range(0, len(a)):
        a1 = a[i][i]
        for ii in range(0, len(a)):
            a2 = a[ii][i]
            if i != ii:
                for iii in range(0, len(a)):
                    a[ii][iii] = a[ii][iii] - a[i][iii]*a2/a1
                    p[ii][iii] = p[ii][iii] - p[i][iii]*a2/a1
    for i in range(0, len(a)):
        a1 = a[i][i]
        for ii in range(0, len(a)):
            p[i][ii] = p[i][ii]/a1
    return p


def main():
    try:
        n = int(input("정방행렬의 차수를 입력하세요: "))
        if n <= 0:
            print("k는 0이상의 정수여야 합니다.")
            return
        print(f"{n}차 정방행렬 A를 입력하세요 (각 행을 공백으로 구분하여 입력하세요):")
        mat = []
        for i in range(n):
            while True:
                try:
                    row_input = input(f"{i+1}번째 행 입력: ").strip()
                    row_values = [float(x) for x in row_input.split()]
                    if len(row_values) != n:
                        print(f"각 행은 {n}개의 숫자로 이루어져야 합니다.")
                        continue
                    mat.append(row_values)
                    break
                except ValueError:
                    print("유효한 숫자로 입력해 주세요.")
        A = mat
        print("입력된 행렬 A:")
        for row in A:
            print(" ".join(f"{val:.2f}" for val in row))
        try:
            inv_A = inverse_with_det(A)
            print("\n행렬식을 이용해서 구한 행렬 A의 역행렬:")
            for row in inv_A:
                print(" ".join(f"{val:.2f}" for val in row))
            inv_A_gj = gauss_jordan_inverse(A)
            print("\n가우스-조던 소거법으로 구한 행렬 A의 역행렬:")
            for row in inv_A_gj:
                print(" ".join(f"{val:.2f}" for val in row))
            if " ".join(f"{val:.2f}" for val in row) == " ".join(f"{val:.2f}" for val in row):
                print("\n두 방법으로 구한 결과가 동일합니다.")
        except has_inverse(A) == False:
            print("이 행렬은 역행렬이 없습니다.")
    except Exception as e:
        print(f"예외 발생: {e}")


# 프로그램 시작
if __name__ == "__main__":
    main()
