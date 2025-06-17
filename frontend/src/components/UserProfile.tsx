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
    fetch(`/api/user/${user_id}`, { credentials: "include" })
      .then(res => res.json())
      .then(data => setUser(data))
      .finally(() => setLoading(false));
  }, [user_id]);

  if (loading) return <div className="text-white p-8">Loading...</div>;
  if (!user || !user.email) return <div className="text-white p-8">User not found.</div>;

  return (
    <div className="max-w-md mx-auto p-8 bg-gray-900 rounded-lg shadow text-white mt-12">
      <h2 className="text-2xl font-bold mb-4">User Profile</h2>
      <div className="mb-2"><b>Email:</b> {user.email}</div>
      <button className="mt-4 px-4 py-2 bg-green-700 rounded" onClick={() => navigate(-1)}>Back</button>
    </div>
  );
};

export default UserProfile; 