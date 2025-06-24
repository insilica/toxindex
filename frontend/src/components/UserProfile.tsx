import React, { useEffect, useState } from "react";
import { useParams, useNavigate } from "react-router-dom";

const UserProfile: React.FC = () => {
  const { user_id } = useParams<{ user_id: string }>();
  const [user, setUser] = useState<any>(null);
  const [loading, setLoading] = useState(true);
  const navigate = useNavigate();

  useEffect(() => {
    if (!user_id) return;
    setLoading(true);
    fetch(`/api/users/${user_id}`, { credentials: "include" })
      .then(res => res.json())
      .then(data => setUser(data))
      .finally(() => setLoading(false));
  }, [user_id]);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!user || !user.email) return <div className="text-white p-8">User not found.</div>;

  return (
    <div className="min-h-screen flex flex-col items-center justify-center bg-gradient-to-br from-gray-950 via-gray-900 to-gray-800 py-16 px-4">
      <div className="w-full max-w-lg bg-gray-900 rounded-2xl shadow-2xl p-10 flex flex-col items-center border border-gray-800">
        <div className="flex flex-col items-center mb-6">
          <div className="w-20 h-20 rounded-full bg-gradient-to-br from-green-400 to-green-700 flex items-center justify-center mb-4 shadow-lg">
            <svg width="48" height="48" fill="none" viewBox="0 0 24 24" stroke="currentColor" strokeWidth="2.2" className="text-white">
              <circle cx="12" cy="8" r="4" />
              <path d="M4 20c0-2.5 3.5-4 8-4s8 1.5 8 4" />
            </svg>
          </div>
          <h2 className="text-3xl font-extrabold text-white mb-1 tracking-tight">User Profile</h2>
          <span className="text-green-400 font-mono text-xs bg-gray-800 rounded px-2 py-1 mt-1 break-all" title={user.user_id}>{user.user_id}</span>
        </div>
        <div className="w-full flex flex-col items-center">
          <div className="mb-2 w-full flex flex-col items-center">
            <span className="uppercase text-xs text-gray-400 tracking-widest mb-1">Email</span>
            <span className="text-lg font-semibold text-white break-all">{user.email}</span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default UserProfile; 