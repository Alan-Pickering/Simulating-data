* Encoding: UTF-8.
INPUT PROGRAM.
+ DO REPEAT sim=1 TO 1000.
+    LOOP subj_id=1 TO 50.
+      COMPUTE simstudy_num=sim.
+      COMPUTE group=1.
+      IF  (subj_id > 25) group=2.
+      COMPUTE DV=RV.NORMAL(0,1).
+      END CASE.
+    END LOOP.
+  END REPEAT.
+  END FILE.
END INPUT PROGRAM.
EXECUTE

*an alternative way to do the same thing.
*incorporating a random number seed
SET RNG=MT MTINDEX=2000000.
INPUT PROGRAM.
+  LOOP #I = 1 TO 200.
+     DO REPEAT simvar = dv1 TO dv1000.
+          COMPUTE group=1.
+           IF  (#I > 100) group=2.
+          COMPUTE subj_id = #I.
+          COMPUTE simvar=RV.NORMAL(0,1).
+     END REPEAT.                    
+     END CASE.                        
+   END LOOP.
+   END FILE.
END INPUT PROGRAM.
EXECUTE

