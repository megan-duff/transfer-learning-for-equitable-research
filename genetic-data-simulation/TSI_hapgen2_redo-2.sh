#!/bin/bash
#SBATCH --job-name=TSI_HM3_GW
#SBATCH --output=TSI_HM3_GW
#SBATCH --error=TSI_HM3_GW_err
#SBATCH -n 16
#SBATCH -p math-alderaan
#SBATCH --array=155,177,265,419,485,550,617,682,847,1061,1210,1211,1280,1281,1282,1283,1284,1285,1288,1289,1290,1291,1292,1293,1294,1295,1296,1297,1298,1299,1300,1301,1302,1303,1304,1305,1306,1307,1308,1309,1311,1312,1313,1314,1315,1316,1317,1318,1319,1320,1321,1322,1323,1327,1328,1329,1330,1331,1332,1333,1334,1335,1336,1337,1338,1339,1340,1341,1342,1343,1344,1345,1346,1347,1348,1349,1350,1351,1353,1354,1355,1356,1357,1358,1359,1360,1361,1362,1363,1364,1365,1366,1367,1368,1369,1370,1371,1372,1374,1375,1376,1377,1378,1379,1380,1381,1382,1383,1384,1385,1386,1387,1389,1390,1391,1392,1393,1394,1395,1396,1397,1398,1399,1400,1401,1402,1403,1404,1405,1407,1409,1410,1411,1412,1413,1414,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1429,1430,1431,1432,1433,1434,1435,1436,1437,1438,1439,1440,1441,1442,1444,1445,1446,1447,1448,1449,1450,1451,1452,1454,1455,1456,1457,1458,1459,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1513,1514,1515,1516,1517,1519,1520,1521,1522,1523,1524,1525,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1547,1548,1549,1550,1551,1552,1553,1554,1555,1556,1557,1558,1559,1560,1561,1562,1563,1564,1565,1566,1567,1568,1569,1570,1571,1572,1573,1574,1575,1576,1577,1578,1579,1580,1581,1582,1583,1584,1585,1586,1587,1588,1589,1590,1591,1592,1593,1594,1595,1596,1597,1598,1599,1600,1601,1602,1603,1604,1605,1606,1607,1608,1609,1610,1611,1612,1613,1614,1615,1616,1617,1618,1619,1620,1621,1622,1623,1624,1625,1626,1627,1628,1629,1630,1631,1632,1633,1634,1635,1636,1637,1638,1639,1640,1641,1642,1643,1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1657,1658,1659,1660,1661,1662,1663,1664,1665,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676,1677,1678,1681,1683,1684,1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702,1703,1704,1705,1706,1707,1709,1710,1711,1712,1713,1714,1715,1716,1717,1718,1719,1720,1721,1722,1723,1724,1725,1726,1728,1729,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747,1748,1749,1750,1751,1752,1753,1754,1755,1756,1757,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1814,1815,1816,1817,1818,1834,1835,1836,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875,1876,1877,1878,1879,1880,1881,1882,1883,1884,1885,1886,1887,1888,1889,1890,1891,1892,1893,1894,1895,1896,1897,1898,1899


#Breaks up 400 simulations of 22 chromosomes into 8400 jobs to be submitted at once

chr=$(expr ${SLURM_ARRAY_TASK_ID} % 22 + 1)
((sim_num=${SLURM_ARRAY_TASK_ID}/22 + 1))

#Create dummy disease locus for HAPGEN2 because even though I am simulating under their control model it won't run wi$

dummyDL1=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr1_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL2=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr2_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL3=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr3_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL4=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr4_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL5=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr5_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL6=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr6_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL7=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr7_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL8=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr8_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL9=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr9_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL10=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr10_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL11=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr11_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL12=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr12_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL13=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr13_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL14=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr14_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL15=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr15_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL16=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr16_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL17=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr17_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL18=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr18_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL19=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr19_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL20=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr20_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL21=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr21_1000G_phase3_HM3.legend|head -40|tail -39)
dummyDL22=$(awk '{print$2}' /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr22_1000G_phase3_HM3.legend|head -40|tail -39)

dummy=dummyDL${chr}
dummy2=(`eval echo $dummy`)
dummy3=(`echo "${!dummy2}"`)

### simulate with HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/hapgen/hapgen2.html for details on HAPGEN2
### see https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html for details on effective sample size es$
### see https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html for details on files used for 1000G Phase 3 data on $

#Simulates 20000 TSI cases for each chromosome

./hapgen2 -m /storage/math/projects/duffme_gensim/Recom_Rates/TSI/TSI-${chr}-final.txt \
        -l /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr${chr}_1000G_phase3_HM3.legend \
        -h /storage/math/projects/duffme_gensim/HAPGEN2_reference_files/TSI/TSI_chr${chr}_1000G_phase3_HM3.haps \
        -o /scratch/duffme_gensim/Simulations/TSI/sim_${sim_num}/TSI_1000GP_hm3_Phase3_chr${chr}_20k \
        -Ne 11418 \
        -n 20000 0 \
	-dl ${dummy3[0]} 1 1 1 ${dummy3[1]} 1 1 1 ${dummy3[2]} 1 1 1 ${dummy3[3]} 1 1 1 ${dummy3[4]} 1 1 1 ${dummy3[5]} 1 1 1 ${dummy3[6]} 1 1 1 ${dummy3[7]} 1 1 1 ${dummy3[8]} 1 1 1 ${dummy3[9]} 1 1 1 ${dummy3[10]} 1 1 1 ${dummy3[11]} 1 1 1 ${dummy3[12]} 1 1 1 ${dummy3[13]} 1 1 1 ${dummy3[14]} 1 1 1 ${dummy3[15]} 1 1 1 ${dummy3[16]} 1 1 1 ${dummy3[17]} 1 1 1 ${dummy3[18]} 1 1 1 ${dummy3[19]} 1 1 1 ${dummy3[20]} 1 1 1 ${dummy3[21]} 1 1 1 ${dummy3[22]} 1 1 1 ${dummy3[23]} 1 1 1 ${dummy3[24]} 1 1 1 ${dummy3[25]} 1 1 1 ${dummy3[26]} 1 1 1  ${dummy3[27]} 1 1 1 ${dummy3[28]} 1 1 1 ${dummy3[29]} 1 1 1 ${dummy3[30]} 1 1 1 ${dummy3[31]} 1 1 1 ${dummy3[32]} 1 1 1 ${dummy3[33]} 1 1 1 ${dummy3[34]} 1 1 1 ${dummy3[35]} 1 1 1 ${dummy3[36]} 1 1 1 ${dummy3[37]} 1 1 1 ${dummy3[38]} 1 1 1 ${dummy3[39]} 1 1 1 \
	-no_haps_output

if ! grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/TSI/sim_${sim_num}/TSI_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
	echo chr: ${chr} sim: ${sim_num} >> TSI_failed_sims.txt
fi

if grep 'Simulating haplotypes ... done' /scratch/duffme_gensim/Simulations/TSI/sim_${sim_num}/TSI_1000GP_hm3_Phase3_chr${chr}_20k.summary
then
cd /scratch/duffme_gensim/Simulations/TSI/sim_${sim_num}
rm *cases*
for chr in {1..22};
do
../../.././plink --data TSI_1000GP_hm3_Phase3_chr${chr}_20k.controls\
		--make-bed\
		--oxford-single-chr ${chr}\
		--allow-no-sex\
		--out TSI_1000GP_hm3_Phase3_chr${chr}_20k.controls;
done
fi
