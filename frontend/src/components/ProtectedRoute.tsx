import React, { useEffect, useState } from "react";
import { Navigate } from "react-router-dom";

const ProtectedRoute: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [auth, setAuth] = useState<null | boolean>(null);

  useEffect(() => {
    fetch("/api/users/me", { credentials: "include", cache: "no-store" })
      .then(res => {
        if (res.status === 200) {
          setAuth(true);
        } else {
          setAuth(false);
        }
      })
      .catch(() => {
        setAuth(false);
      });
  }, []);

  if (auth === null) {
    return <div>Loading...</div>;
  }

  return auth ? <>{children}</> : <Navigate to="/login" />;
};

export default ProtectedRoute; 