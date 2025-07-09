import React from 'react';
import { useAdmin } from '../../hooks/useAdmin';

interface AdminRouteProps {
  children: React.ReactNode;
  fallback?: React.ReactNode;
}

const AdminRoute: React.FC<AdminRouteProps> = ({ 
  children, 
  fallback = (
    <div className="flex justify-center items-center h-64">
      <div className="text-white text-center">
        <div className="text-xl font-bold mb-2">Access Denied</div>
        <div className="text-gray-300">You need administrator privileges to access this page.</div>
      </div>
    </div>
  ) 
}) => {
  const { isAdmin, loading } = useAdmin();

  if (loading) {
    return (
      <div className="flex justify-center items-center h-64">
        <div className="text-white">Loading...</div>
      </div>
    );
  }

  if (!isAdmin) {
    return <>{fallback}</>;
  }

  return <>{children}</>;
};

export default AdminRoute; 