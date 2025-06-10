import React, { useEffect, useState } from "react";
import { useNavigate } from "react-router-dom";

const Layout: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const navigate = useNavigate();
  const [auth, setAuth] = useState<null | boolean>(null);

  useEffect(() => {
    fetch("/api/me", { credentials: "include", cache: "no-store" })
      .then(res => {
        if (res.status === 200) setAuth(true);
        else setAuth(false);
      })
      .catch(() => setAuth(false));
  }, []);

  const handleLogout = () => {
    fetch("/api/logout", {
      method: "GET",
      credentials: "include",
    }).then(() => {
      navigate("/login");
    });
  };

  return (
    <div className="flex min-h-screen w-screen" style={{ background: '#10291a' }}>
      {/* Sidebar */}
      <aside className="w-64 bg-black border-r shadow-sm p-4 flex flex-col">
        <h1 className="text-2xl font-bold mb-6 text-white">Your App</h1>
        {/* Add sidebar navigation here */}
        <nav>
          <ul>
            <li className="mb-2"><a href="#" className="hover:text-blue-600 text-white">Environments</a></li>
            <li className="mb-2"><a href="#" className="hover:text-blue-600 text-white">Tasks</a></li>
            {/* Add more links as needed */}
          </ul>
        </nav>
        {/* Only show logout if authenticated */}
        {auth && (
          <button
            onClick={handleLogout}
            className="mt-8 px-4 py-2 bg-red-500 text-white rounded hover:bg-red-600 w-full"
          >
            Logout
          </button>
        )}
      </aside>
      {/* Main Content */}
      <main className="flex-1 flex items-center justify-center min-h-screen p-8">{children}</main>
    </div>
  );
};

export default Layout; 