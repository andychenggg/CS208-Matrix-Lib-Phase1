//#include <stdio.h>
//#include "matrix.h"
//#include "string.h"
//int main()
//{
//
//    float f[8] = {1.f,0.f,0.f,0.f,
//                  0.f,2.f,0.f,0.f};
//    Matrix *p = createMatrix_with_pFloat(2,4,f,8);
//    determinant(p);
//
//    char c [4];
//    strlen(c)
//}

//#include <stdio.h>
//#include <stdlib.h>
//#include <libpq-fe.h>
//
//int main() {
//    PGconn *conn = PQconnectdb("host=localhost port=5432 dbname=mydb user=myuser password=mypass");
//    if (PQstatus(conn) != CONNECTION_OK) {
//        fprintf(stderr, "Connection to database failed: %s", PQerrorMessage(conn));
//        PQfinish(conn);
//        exit(1);
//    }
//
//    PGresult *res = PQexec(conn, "SELECT * FROM mytable");
//    if (PQresultStatus(res) != PGRES_TUPLES_OK) {
//        fprintf(stderr, "Query failed: %s", PQerrorMessage(conn));
//        PQclear(res);
//        PQfinish(conn);
//        exit(1);
//    }
//
//    int rows = PQntuples(res);
//    int cols = PQnfields(res);
//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < cols; j++) {
//            printf("%s\t", PQgetvalue(res, i, j));
//        }
//        printf("\n");
//    }
//
//    PQclear(res);
//    PQfinish(conn);
//    return 0;
//}
